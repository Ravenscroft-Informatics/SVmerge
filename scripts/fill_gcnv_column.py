#!/usr/bin/env python3
"""
Convenience script tying Snakefile.gcnv's output back into the merge
pipeline: fills in (or adds) the "gcnv_vcf" column of the merge pipeline's
samples.tsv from Snakefile.gcnv's per-sample genotyped-segments VCFs.

Not part of either Snakemake DAG on purpose -- this is a one-shot,
run-it-yourself step between "the gCNV workflow finished" and "point
config.yaml's enable_gcnv at samples.tsv", not something that needs to
re-run automatically every time either pipeline executes.

Matches samples by sample_id: any merge-pipeline sample without a
corresponding {sample_id}.genotyped-segments.vcf.gz under gcnv_outdir is
left with an empty gcnv_vcf value and reported at the end, rather than
silently dropped or defaulted -- if you're running enable_gcnv, every
sample needs one.

Usage:
    python3 fill_gcnv_column.py <merge_samples_tsv> <gcnv_outdir> <output_tsv>

<gcnv_outdir> is Snakefile.gcnv's config_gcnv.yaml "outdir" value -- this
script looks for <gcnv_outdir>/genotyped/<sample_id>.genotyped-segments.vcf.gz
for each sample.
"""
import os
import sys

import pandas as pd


def main():
    if len(sys.argv) != 4:
        sys.exit(f"usage: {sys.argv[0]} <merge_samples_tsv> <gcnv_outdir> <output_tsv>")
    merge_samples_tsv, gcnv_outdir, out_path = sys.argv[1:4]

    df = pd.read_csv(merge_samples_tsv, sep="\t", dtype=str)
    if "sample_id" not in df.columns:
        sys.exit(f"ERROR: {merge_samples_tsv} has no 'sample_id' column.")

    genotyped_dir = os.path.join(gcnv_outdir, "genotyped")
    missing = []

    def _lookup(sample_id):
        vcf_path = os.path.join(genotyped_dir, f"{sample_id}.genotyped-segments.vcf.gz")
        if os.path.isfile(vcf_path):
            return vcf_path
        missing.append(sample_id)
        return ""

    df["gcnv_vcf"] = df["sample_id"].map(_lookup)

    df.to_csv(out_path, sep="\t", index=False)

    n_total = len(df)
    n_found = n_total - len(missing)
    print(f"[fill_gcnv_column] {n_found}/{n_total} samples matched a gCNV VCF.", file=sys.stderr)
    if missing:
        print(
            f"[fill_gcnv_column] WARNING: no genotyped-segments.vcf.gz found for "
            f"{len(missing)} sample(s), left blank in gcnv_vcf: {missing}",
            file=sys.stderr,
        )
        print(
            "[fill_gcnv_column] Check these ran through Snakefile.gcnv (or were "
            "included in its samples_gcnv.csv) before enabling enable_gcnv on them.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()