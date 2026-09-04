#!/usr/bin/env python3
"""
Build the sample ploidy table GATK SVCluster requires (--ploidy-table):
tab-delimited, header "SAMPLE" + one column per contig, each row a sample's
integer ploidy per contig.

Autosomes are always ploidy 2. Sex chromosomes are derived from a required
"sex" column in samples.tsv (M/F, case-insensitive; "male"/"female" also
accepted): female -> X=2, Y=0; male -> X=1, Y=1.

This is the simple, standard-human-germline-WGS case. It will NOT match a
gCNV run that used anything more specific (e.g. mosaicism handling,
non-standard ploidy regions, non-human samples) -- if DetermineGermlineContigPloidy
was run with custom ploidy priors, prefer reusing its own output over this
script, since the SVCluster ploidy table should match what the gCNV calls
were actually made under.

Usage:
    python3 make_ploidy_table.py <samples_tsv> <comma_separated_contigs> <output_tsv>
"""
import sys

import pandas as pd


def classify_contig(contig: str) -> str:
    c = contig.lower()
    if c in ("chrx", "x"):
        return "X"
    if c in ("chry", "y"):
        return "Y"
    return "autosome"


def main():
    if len(sys.argv) != 4:
        sys.exit(f"usage: {sys.argv[0]} <samples_tsv> <comma_separated_contigs> <output_tsv>")
    samples_tsv, contigs_arg, out_path = sys.argv[1:4]
    contigs = [c for c in contigs_arg.split(",") if c]
    if not contigs:
        sys.exit("ERROR: no contigs given.")

    # encoding="utf-8-sig" strips a leading UTF-8 BOM if present (harmless if
    # not); column names are stripped of surrounding whitespace too, since a
    # stray trailing \r or space on the header line would otherwise make a
    # column look present-but-not-matching -- see make_contig_ploidy_priors.py
    # for the same hardening, added there after hitting exactly this.
    df = pd.read_csv(samples_tsv, sep="\t", dtype=str, encoding="utf-8-sig")
    df.columns = [c.strip() for c in df.columns]
    if "sample_id" not in df.columns:
        sys.exit(f"ERROR: samples.tsv has no 'sample_id' column. Columns found: {list(df.columns)!r}")
    if "sex" not in df.columns:
        sys.exit(
            "ERROR: samples.tsv has no 'sex' column. make_ploidy_table.py needs "
            "one M/F value per sample to assign X/Y ploidy. Add a 'sex' column "
            f"(M or F, case-insensitive) to samples.tsv. Columns found: {list(df.columns)!r}"
        )

    rows = []
    bad = []
    for _, row in df.iterrows():
        sample = row["sample_id"]
        sex_raw = str(row["sex"]).strip().upper()
        if sex_raw in ("M", "MALE"):
            sex = "M"
        elif sex_raw in ("F", "FEMALE"):
            sex = "F"
        else:
            bad.append((sample, row["sex"]))
            continue

        ploidy = {}
        for contig in contigs:
            kind = classify_contig(contig)
            if kind == "autosome":
                ploidy[contig] = 2
            elif kind == "X":
                ploidy[contig] = 2 if sex == "F" else 1
            else:  # "Y"
                ploidy[contig] = 0 if sex == "F" else 1
        rows.append((sample, ploidy))

    if bad:
        detail = "; ".join(f"{s}='{v}'" for s, v in bad)
        sys.exit(
            f"ERROR: unrecognized 'sex' value(s) for {len(bad)} sample(s): {detail}. "
            "Expected M/F (or male/female), case-insensitive."
        )

    with open(out_path, "w") as fh:
        fh.write("SAMPLE\t" + "\t".join(contigs) + "\n")
        for sample, ploidy in rows:
            fh.write(sample + "\t" + "\t".join(str(ploidy[c]) for c in contigs) + "\n")

    print(
        f"[make_ploidy_table] wrote ploidy table for {len(rows)} samples, "
        f"{len(contigs)} contigs, to {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()