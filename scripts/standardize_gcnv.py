#!/usr/bin/env python3
"""
Standardizes raw GATK PostprocessGermlineCNVCalls (gCNV) VCF output for
GATK SVCluster:

  - INFO/ALGORITHMS = "depth", if not already present (gCNV is purely a
    depth-based caller; SVCluster uses ALGORITHMS to distinguish
    depth-only/PESR/mixed comparisons and pick the right tolerance set).
  - FORMAT/ECN (expected copy number, i.e. ploidy), looked up per
    sample+contig from a ploidy table (see make_ploidy_table.py), if not
    already present.
  - INFO/SVTYPE, derived from FORMAT/CN vs FORMAT/ECN (CN < ECN -> DEL,
    CN > ECN -> DUP) *only* if not already present in the record --
    raw gCNV segments VCFs typically already carry SVTYPE, so this is a
    defensive fallback rather than the expected path. Records where
    CN == ECN (no copy-number change) and SVTYPE can't otherwise be
    determined are left without SVTYPE; downstream filtering/SVCluster
    can drop them.

Assumes single-sample input, matching PostprocessGermlineCNVCalls' own
per-sample VCF output.

Usage:
    python3 standardize_gcnv.py <input.vcf[.gz]> <output.vcf> <ploidy_table.tsv> <sample_id>
"""
import sys

import pysam


def load_ploidy(ploidy_table_path, sample_id):
    with open(ploidy_table_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if not header or header[0] != "SAMPLE":
            sys.exit(
                f"ERROR: {ploidy_table_path} does not look like a ploidy table "
                "(expected 'SAMPLE' as the first header column)."
            )
        contigs = header[1:]
        for line in fh:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == sample_id:
                return dict(zip(contigs, (int(x) for x in fields[1:])))
    sys.exit(f"ERROR: sample '{sample_id}' not found in ploidy table {ploidy_table_path}.")


def main():
    if len(sys.argv) != 5:
        sys.exit(
            f"usage: {sys.argv[0]} <input.vcf[.gz]> <output.vcf> "
            "<ploidy_table.tsv> <sample_id>"
        )
    in_path, out_path, ploidy_table_path, sample_id = sys.argv[1:5]

    ploidy = load_ploidy(ploidy_table_path, sample_id)

    vin = pysam.VariantFile(in_path)
    if sample_id not in vin.header.samples:
        sys.exit(
            f"ERROR: sample '{sample_id}' not found in {in_path}'s header "
            f"(found: {list(vin.header.samples)}). standardize_gcnv.py expects "
            "single-sample PostprocessGermlineCNVCalls output."
        )

    header = pysam.VariantHeader()
    for line in str(vin.header).rstrip("\n").split("\n"):
        if not line or line.startswith("#CHROM"):
            continue
        header.add_line(line)
    for s in vin.header.samples:
        header.add_sample(s)

    if "ALGORITHMS" not in header.info:
        header.info.add(
            "ALGORITHMS", ".", "String", "Source algorithms/callers supporting this call"
        )
    if "ECN" not in header.formats:
        header.formats.add("ECN", 1, "Integer", "Expected copy number (ploidy)")
    if "SVTYPE" not in header.info:
        header.info.add("SVTYPE", 1, "String", "Type of structural variant")

    vout = pysam.VariantFile(out_path, "w", header=header)

    n_algorithms_added = n_ecn_added = n_svtype_derived = n_svtype_undetermined = 0

    for rec in vin:
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
            try:
                new_rec.info[key] = value
            except (TypeError, OverflowError):
                pass

        if "ALGORITHMS" not in rec.info:
            new_rec.info["ALGORITHMS"] = "depth"
            n_algorithms_added += 1

        ecn = ploidy.get(rec.contig)
        gt_fields = rec.samples[sample_id]
        cn = None
        if "CN" in gt_fields.keys():
            cn_val = gt_fields["CN"]
            cn = cn_val[0] if isinstance(cn_val, tuple) else cn_val

        svtype = new_rec.info.get("SVTYPE")
        if svtype is None:
            if cn is not None and ecn is not None and cn != ecn:
                new_rec.info["SVTYPE"] = "DEL" if cn < ecn else "DUP"
                n_svtype_derived += 1
            else:
                n_svtype_undetermined += 1

        for key in gt_fields.keys():
            try:
                new_rec.samples[sample_id][key] = gt_fields[key]
            except TypeError:
                pass
        if ("ECN" not in gt_fields.keys() or gt_fields.get("ECN") is None) and ecn is not None:
            new_rec.samples[sample_id]["ECN"] = ecn
            n_ecn_added += 1

        vout.write(new_rec)

    vin.close()
    vout.close()

    print(
        f"[standardize_gcnv] {sample_id}: ALGORITHMS added={n_algorithms_added}, "
        f"ECN added={n_ecn_added}, SVTYPE derived={n_svtype_derived}, "
        f"SVTYPE undetermined (CN==ECN or missing CN/ploidy)={n_svtype_undetermined}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()