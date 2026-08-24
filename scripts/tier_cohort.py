#!/usr/bin/env python3
"""
Stage 2 post-processing: after `svdb --merge` combines all per-sample
merged VCFs into one cohort VCF, this script:

  1. Drops SVDB's verbose per-input-file provenance INFO fields (the
     "<filename>_INFO", "<filename>_SAMPLE", ... dump SVDB adds per
     contributing input file).
  2. Aggregates the per-sample FORMAT/CS tag (added at stage 1 by
     add_caller_support.py) across all samples at each site, and adds
     cohort-level INFO fields summarizing caller support:
       N_CALLED       - samples with a non-missing genotype at this site
       N_BOTH         - of those, called by both Manta and Dysgu
       N_MANTA_ONLY   - of those, called by Manta only
       N_DYSGU_ONLY   - of those, called by Dysgu only
       CALLER_TIER    - DUAL_CALLER_MAJORITY / MANTA_MAJORITY /
                        DYSGU_MAJORITY / MIXED / NO_CALLS

Usage:
    python3 tier_cohort.py input.vcf output.vcf
"""
import re
import sys
from collections import Counter

import pysam

# Matches SVDB's per-input-file provenance fields, e.g.
# "sample1.merged.sorted_INFO", "sample42_SAMPLE", "s3.tagged_CHROM", ...
DROP_INFO_PATTERN = re.compile(
    r"^.+_(INFO|SAMPLE|CHROM|POS|QUAL|FILTERS)$"
)
# Also drop these SVDB bookkeeping fields if present (harmless to keep,
# but redundant with the counts we compute below).
DROP_INFO_EXACT = {"VARID", "set", "FOUNDBY", "svdb_origin"}


def classify(n_both: int, n_manta: int, n_dysgu: int) -> str:
    n_called = n_both + n_manta + n_dysgu
    if n_called == 0:
        return "NO_CALLS"
    if n_both / n_called >= 0.5:
        return "DUAL_CALLER_MAJORITY"
    if n_manta > n_dysgu:
        return "MANTA_MAJORITY"
    if n_dysgu > n_manta:
        return "DYSGU_MAJORITY"
    return "MIXED"


def main(in_path: str, out_path: str) -> None:
    vin = pysam.VariantFile(in_path)
    samples = list(vin.header.samples)

    info_ids_to_drop = {
        key
        for key in list(vin.header.info)
        if DROP_INFO_PATTERN.match(key) or key in DROP_INFO_EXACT
    }
    drop_header_line = (
        re.compile(
            r"^##INFO=<ID=("
            + "|".join(re.escape(k) for k in info_ids_to_drop)
            + r"),"
        )
        if info_ids_to_drop
        else None
    )

    header = pysam.VariantHeader()
    for line in str(vin.header).rstrip("\n").split("\n"):
        if not line or line.startswith("#CHROM"):
            continue
        if drop_header_line and drop_header_line.match(line):
            continue
        header.add_line(line)
    for s in samples:
        header.add_sample(s)

    header.info.add("N_CALLED", 1, "Integer", "Samples with a non-missing genotype at this site")
    header.info.add("N_BOTH", 1, "Integer", "Called samples supported by both Manta and Dysgu")
    header.info.add("N_MANTA_ONLY", 1, "Integer", "Called samples supported by Manta only")
    header.info.add("N_DYSGU_ONLY", 1, "Integer", "Called samples supported by Dysgu only")
    header.info.add(
        "CALLER_TIER",
        1,
        "String",
        "Site-level caller-support tier: DUAL_CALLER_MAJORITY, MANTA_MAJORITY, "
        "DYSGU_MAJORITY, MIXED, or NO_CALLS",
    )

    vout = pysam.VariantFile(out_path, "w", header=header)
    tier_counts = Counter()

    for rec in vin:
        cs_counts = Counter()
        for s in samples:
            gt = rec.samples[s].get("GT")
            is_called = gt is not None and any(a is not None for a in gt) and gt != (0, 0)
            if not is_called:
                continue
            cs = rec.samples[s].get("CS")
            cs_counts[cs] += 1

        n_both = cs_counts.get("BOTH", 0)
        n_manta = cs_counts.get("MANTA_ONLY", 0)
        n_dysgu = cs_counts.get("DYSGU_ONLY", 0)
        n_called = n_both + n_manta + n_dysgu
        tier = classify(n_both, n_manta, n_dysgu)
        tier_counts[tier] += 1

        new_rec = vout.new_record(
            contig=rec.contig,
            start=rec.start,
            stop=rec.stop,
            alleles=rec.alleles,
            id=rec.id,
            qual=rec.qual,
            filter=list(rec.filter.keys()) if rec.filter else None,
        )
        for key, value in rec.info.items():
            if key in info_ids_to_drop:
                continue
            new_rec.info[key] = value
        new_rec.info["N_CALLED"] = n_called
        new_rec.info["N_BOTH"] = n_both
        new_rec.info["N_MANTA_ONLY"] = n_manta
        new_rec.info["N_DYSGU_ONLY"] = n_dysgu
        new_rec.info["CALLER_TIER"] = tier

        for s in samples:
            for key in rec.samples[s].keys():
                new_rec.samples[s][key] = rec.samples[s][key]

        vout.write(new_rec)

    vin.close()
    vout.close()

    total = sum(tier_counts.values())
    summary = ", ".join(f"{k}={v}" for k, v in sorted(tier_counts.items()))
    print(f"[tier_cohort] {total} sites: {summary}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <input.vcf[.gz]> <output.vcf>")
    main(sys.argv[1], sys.argv[2])
