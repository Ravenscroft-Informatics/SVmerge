#!/usr/bin/env python3
"""
Stage 1 post-processing: after `svdb --merge` combines a single sample's
Manta + Dysgu calls, this script:

  1. Reads the SUPP_VEC INFO field SVDB adds (a binary string like "11",
     "10", "01" indicating which of the two input files support each
     merged record; input file order matches the --vcf order given to
     svdb --merge).
  2. Translates SUPP_VEC into a per-sample FORMAT/CS tag:
       11 -> BOTH
       10 -> MANTA_ONLY   (assumes manta was --vcf file 1)
       01 -> DYSGU_ONLY   (assumes dysgu was --vcf file 2)
  3. Drops SVDB's verbose per-file provenance INFO fields (the
     "1_INFO", "2_SAMPLE", "1_CHROM", ... dump, plus VARID/set/FOUNDBY/
     svdb_origin/SUPP_VEC) so they don't collide with the INFO fields
     SVDB adds again at the cohort-merge stage.

Usage:
    python3 add_caller_support.py input.vcf output.vcf

Assumes the input VCF is single-sample and was produced with:
    svdb --merge --vcf manta.vcf.gz:1 dysgu.vcf.gz:2 --priority 1,2 ...
so that SUPP_VEC position 1 = manta, position 2 = dysgu.
"""
import re
import sys

import pysam

SUPP_VEC_TO_CS = {
    "11": "BOTH",
    "10": "MANTA_ONLY",
    "01": "DYSGU_ONLY",
}

# INFO fields SVDB adds that we don't want to carry into the cohort merge
# (they'd collide with the new fields SVDB adds at that stage, since the
# field names are reused). Matches things like 1_INFO, 2_SAMPLE, 1_CHROM,
# VARID, set, FOUNDBY, svdb_origin, SUPP_VEC.
DROP_INFO_PATTERN = re.compile(
    r"^(\d+_(INFO|SAMPLE|CHROM|POS|QUAL|FILTERS)|VARID|set|FOUNDBY|svdb_origin|SUPP_VEC)$"
)


def main(in_path: str, out_path: str) -> None:
    vin = pysam.VariantFile(in_path)

    if len(vin.header.samples) != 1:
        sys.exit(
            f"error: expected a single-sample VCF, found "
            f"{len(vin.header.samples)} samples in {in_path}"
        )
    sample = vin.header.samples[0]

    # Rebuild the header line-by-line rather than copying, so the noisy
    # provenance INFO fields we're about to drop from every record don't
    # linger as unused (and warning-triggering) header declarations.
    info_ids_to_drop = [
        key for key in list(vin.header.info) if DROP_INFO_PATTERN.match(key)
    ]
    drop_header_line = re.compile(
        r'^##INFO=<ID=(' + "|".join(re.escape(k) for k in info_ids_to_drop) + r'),'
    ) if info_ids_to_drop else None

    header = pysam.VariantHeader()
    for line in str(vin.header).rstrip("\n").split("\n"):
        if not line or line.startswith("#CHROM"):
            continue
        if drop_header_line and drop_header_line.match(line):
            continue
        header.add_line(line)
    header.formats.add(
        "CS",
        1,
        "String",
        "Caller support for this sample's genotype: BOTH, MANTA_ONLY, "
        "DYSGU_ONLY, or UNKNOWN if SUPP_VEC was missing/unrecognized.",
    )
    for s in vin.header.samples:
        header.add_sample(s)

    vout = pysam.VariantFile(out_path, "w", header=header)

    n_both = n_manta = n_dysgu = n_unknown = 0
    n_format_dropped = {}  # key -> count of records where copy failed


    for rec in vin:
        supp_vec = rec.info.get("SUPP_VEC")
        cs = SUPP_VEC_TO_CS.get(supp_vec, "UNKNOWN")
        if cs == "BOTH":
            n_both += 1
        elif cs == "MANTA_ONLY":
            n_manta += 1
        elif cs == "DYSGU_ONLY":
            n_dysgu += 1
        else:
            n_unknown += 1

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
        new_rec.samples[sample]["GT"] = rec.samples[sample]["GT"]
        new_rec.samples[sample].phased = rec.samples[sample].phased
        for key in rec.samples[sample].keys():
            if key == "GT":
                continue
            try:
                new_rec.samples[sample][key] = rec.samples[sample][key]
            except TypeError:
                # SVDB sometimes merges a FORMAT field present in both source
                # VCFs under the same name but writes more values into the
                # record than the (shared) header Number declaration allows
                # -- e.g. a nominally Number=1 field ends up encoded with 2
                # values. Drop just this field for this record rather than
                # crash; GT and CS (what this script actually produces) are
                # unaffected.
                n_format_dropped[key] = n_format_dropped.get(key, 0) + 1
        new_rec.samples[sample]["CS"] = cs

        vout.write(new_rec)

    vin.close()
    vout.close()

    total = n_both + n_manta + n_dysgu + n_unknown
    print(
        f"[add_caller_support] {sample}: {total} records "
        f"(BOTH={n_both}, MANTA_ONLY={n_manta}, DYSGU_ONLY={n_dysgu}, "
        f"UNKNOWN={n_unknown})",
        file=sys.stderr,
    )
    if n_format_dropped:
        print(
            "[add_caller_support] dropped malformed FORMAT values "
            "(header/data Number mismatch) in "
            + ", ".join(f"{k}={v} records" for k, v in n_format_dropped.items()),
            file=sys.stderr,
        )

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <input.vcf[.gz]> <output.vcf>")
    main(sys.argv[1], sys.argv[2])
