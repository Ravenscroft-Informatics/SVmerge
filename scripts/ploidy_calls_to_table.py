#!/usr/bin/env python3
"""
Converts DetermineGermlineContigPloidy's own `ploidy-calls` output directory
into the flat --ploidy-table GATK SVCluster expects, as an alternative to
scripts/make_ploidy_table.py's samples.tsv "sex"-column-derived table.

Why prefer this over make_ploidy_table.py, when it's available: make_ploidy_table.py
*assumes* ploidy from samples.tsv's "sex" column (autosome always 2; X/Y from a
simple M/F rule) -- it has no way to know about real sex-chromosome
aneuploidies (XXY, XO, mosaic loss, etc.) or a wrong "sex" value in samples.tsv.
DetermineGermlineContigPloidy's own calls are the *actual* per-contig ploidy
the tool inferred from each sample's real coverage data, and are exactly what
GermlineCNVCaller/PostprocessGermlineCNVCalls used as the baseline copy number
when generating the gCNV VCFs -- reusing them here guarantees SVCluster's
ploidy table agrees with what the gCNV calls were actually made under. This
is GATK's own recommendation (see make_ploidy_table.py's docstring, which
already flags this trade-off).

Input layout -- Snakefile.gcnv's determine_germline_contig_ploidy rule
produces this at `{gcnv outdir}/ploidy/ploidy-calls/`:

    ploidy-calls/
      SAMPLE_0/
        contig_ploidy.tsv      <- "@RG\tID:GATKCopyNumber\tSM:<sample_name>"
                                   comment line, then a CONTIG/PLOIDY/PLOIDY_GQ
                                   header, then one data row per contig
        global_read_depth.tsv
        sample_name.txt
        ...
      SAMPLE_1/
        ...

(format confirmed against gcnvkernel's own writer: io_ploidy.py's
SamplePloidyWriter._write_sample_contig_ploidy + io_consts.py's column-name
constants, both in broadinstitute/gatk on GitHub.)

SAMPLE_<i> directories are numbered in DetermineGermlineContigPloidy's own
cohort order (0-based) -- in this project, that's Snakefile.gcnv's SAMPLES
list order (samples_gcnv.csv row order). This script does NOT rely on that
index alone lining up with anything: it reads each sample's actual name from
contig_ploidy.tsv's own "@RG...SM:<name>" comment line (falling back to the
sibling sample_name.txt if that line is somehow missing), so a re-ordered or
edited samples_gcnv.csv/samples.tsv can't silently mislabel a sample here.

Output: a single tab-delimited file, header "SAMPLE" + one column per
requested contig, one row per sample, integer PLOIDY values -- the exact
format GATK SVCluster's --ploidy-table requires (see make_ploidy_table.py's
docstring for the same spec). PLOIDY_GQ (genotype quality) is intentionally
dropped: SVCluster's ploidy table has no column for it, and
DetermineGermlineContigPloidy already applied its own --ploidy-gq-filter-min
threshold before deciding each call, so there's nothing left to filter on
here.

Usage:
    python3 ploidy_calls_to_table.py <ploidy_calls_dir> <comma_separated_contigs> <output_tsv>

<ploidy_calls_dir> is Snakefile.gcnv's `{outdir}/ploidy/ploidy-calls` directory
(config_gcnv.yaml's "outdir" + "/ploidy/ploidy-calls").

No pandas dependency on purpose -- this can run in envs/svmerge.yml (or any
plain-Python env) without dragging in the gatk4/pandas version-pin issue
make_contig_ploidy_priors.py's env split (envs/python3.yml) was created to
route around.
"""
import glob
import os
import re
import sys

SAMPLE_NAME_RE = re.compile(r"^@RG.*SM:([^\t\r\n]*)")


def _read_sample_name(sample_dir, contig_ploidy_lines):
    for line in contig_ploidy_lines:
        if line.startswith("@"):
            m = SAMPLE_NAME_RE.match(line.strip())
            if m and m.group(1):
                return m.group(1)
    # Fall back to sample_name.txt, written alongside contig_ploidy.tsv by
    # the same DetermineGermlineContigPloidy run -- only reached if the
    # @RG...SM: comment line is missing or unparseable.
    name_file = os.path.join(sample_dir, "sample_name.txt")
    if os.path.isfile(name_file):
        with open(name_file) as fh:
            name = fh.read().strip()
            if name:
                return name
    return None


def _read_contig_ploidy(sample_dir):
    tsv_path = os.path.join(sample_dir, "contig_ploidy.tsv")
    if not os.path.isfile(tsv_path):
        sys.exit(
            f"ERROR: {tsv_path} not found -- is {sample_dir} really a "
            "DetermineGermlineContigPloidy ploidy-calls/SAMPLE_<i> directory?"
        )
    with open(tsv_path) as fh:
        lines = fh.readlines()

    sample_name = _read_sample_name(sample_dir, lines)
    if sample_name is None:
        sys.exit(
            f"ERROR: could not determine sample name for {sample_dir} "
            "(no '@RG...SM:' comment line in contig_ploidy.tsv, and no "
            "usable sample_name.txt alongside it)."
        )

    data_lines = [line for line in lines if not line.startswith("@")]
    if not data_lines:
        sys.exit(f"ERROR: {tsv_path} has no data rows.")
    header = data_lines[0].rstrip("\n").split("\t")
    try:
        contig_idx = header.index("CONTIG")
        ploidy_idx = header.index("PLOIDY")
    except ValueError:
        sys.exit(
            f"ERROR: {tsv_path}'s header {header!r} doesn't look like a "
            "contig_ploidy.tsv (expected CONTIG and PLOIDY columns)."
        )

    ploidy_by_contig = {}
    for line in data_lines[1:]:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        # DetermineGermlineContigPloidy writes ploidy as an integer already,
        # but round-then-int defensively in case a different gcnvkernel
        # version ever writes it as e.g. "2.0".
        ploidy_by_contig[fields[contig_idx]] = int(round(float(fields[ploidy_idx])))

    return sample_name, ploidy_by_contig


def main():
    if len(sys.argv) != 4:
        sys.exit(f"usage: {sys.argv[0]} <ploidy_calls_dir> <comma_separated_contigs> <output_tsv>")
    ploidy_calls_dir, contigs_arg, out_path = sys.argv[1:4]
    contigs = [c for c in contigs_arg.split(",") if c]
    if not contigs:
        sys.exit("ERROR: no contigs given.")

    sample_dirs = sorted(
        (d for d in glob.glob(os.path.join(ploidy_calls_dir, "SAMPLE_*")) if os.path.isdir(d)),
        key=lambda d: int(os.path.basename(d).split("_", 1)[1]),
    )
    if not sample_dirs:
        sys.exit(
            f"ERROR: no SAMPLE_* subdirectories found under {ploidy_calls_dir} -- "
            "expected Snakefile.gcnv's '{outdir}/ploidy/ploidy-calls' directory."
        )

    rows = []
    seen_names = set()
    for sample_dir in sample_dirs:
        sample_name, ploidy_by_contig = _read_contig_ploidy(sample_dir)
        if sample_name in seen_names:
            sys.exit(
                f"ERROR: duplicate sample name '{sample_name}' encountered "
                f"(from {sample_dir}) -- refusing to silently overwrite a row."
            )
        seen_names.add(sample_name)

        missing_contigs = [c for c in contigs if c not in ploidy_by_contig]
        if missing_contigs:
            sys.exit(
                f"ERROR: {sample_dir}'s contig_ploidy.tsv has no ploidy call for "
                f"{missing_contigs}. These contigs must have been part of the "
                "--contig-ploidy-priors table (and interval list) "
                "DetermineGermlineContigPloidy was run with -- check that the "
                "contigs you're passing here match config_gcnv.yaml's 'contigs' "
                "list used for that run."
            )
        rows.append((sample_name, ploidy_by_contig))

    with open(out_path, "w") as fh:
        fh.write("SAMPLE\t" + "\t".join(contigs) + "\n")
        for sample_name, ploidy_by_contig in rows:
            fh.write(sample_name + "\t" + "\t".join(str(ploidy_by_contig[c]) for c in contigs) + "\n")

    print(
        f"[ploidy_calls_to_table] wrote ploidy table for {len(rows)} samples, "
        f"{len(contigs)} contigs, to {out_path} (source: {ploidy_calls_dir})",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()