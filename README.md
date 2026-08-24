# SVmerge

A Snakemake workflow implementing the two-stage merge framework: per-sample Manta+Dysgu merge, then cohort-level cross-sample merge, with caller-support tiering carried through both stages. Manta's raw output is first passed through an inversion-representation fix and an SV-representation standardization pass before any merging happens. Built and tested against synthetic VCFs and real Manta/Dysgu example records; see "What was validated" below.

## What this does

```
Manta only:
  convertInversion.py (local copy of Manta's own script) converts
  Manta's BND-pair representation of inversions into SVTYPE=INV records

For each of N samples (parallel):
  normalize Manta + Dysgu VCFs
    -> standardize_sv.py: harmonize REF/ALT/SVLEN/END representation across callers, drop BND/TRA records with a decoy/unplaced mate
    -> svdb merge (data_profile=sample: strict, same-individual matching)
    -> tag each record with FORMAT/CS (BOTH / MANTA_ONLY / DYSGU_ONLY)

Scattered by contig (parallel):
  subset every sample's merged VCF to the contig
    -> svdb merge across all samples (data_profile=cohort: permissive, cross-individual matching)
    -> aggregate FORMAT/CS across samples into site-level INFO/CALLER_TIER

Then:
  concat contigs -> annotate AF / F_MISSING
    -> filter (FILTER=PASS, F_MISSING, AF) -> cohort.final.vcf.gz
```

## Requirements

```bash
pip install snakemake pysam pandas svdb
# bcftools + tabix (htslib) — via conda/mamba or your package manager
```

Everything here was tested with `svdb==2.12.0`, `bcftools==1.19`,
`snakemake` (8.x), `pysam==0.24.0`.

Two conda envs are used, activated per-rule via Snakemake's `conda:`
directive (run with `snakemake --use-conda ...`):

- `envs/svmerge.yml` — Python 3 + pysam + bcftools/tabix. Used by essentially every rule (`standardize_sv.py`, `add_caller_support.py`, `tier_cohort.py`, all the `bcftools`/`svdb` steps).
- `envs/manta_convert.yml` — Python 2.7 + samtools + htslib, used *only* by `convert_manta_inversions`. Kept separate because Manta's own `convertInversion.py` is Python 2 only (confirmed by reading the
  script: it does byte/str operations that raise `TypeError` under  Python 3), and conda envs pin a single Python major version.

## Setup

1. Fill in `samples.tsv` (see `samples.tsv.example`):

   ```
   sample_id	manta_vcf	dysgu_vcf
   sample001	/path/to/sample001.manta.vcf.gz	/path/to/sample001.dysgu.vcf.gz
   ...
   ```

2. Edit `config.yaml`:
   - `reference`: path to the reference FASTA used for both callers (needs a `.fai` alongside it — `samtools faidx reference.fasta`).
   - `symbolic_threshold`: size (bp) at or above which `standardize_sv.py` rewrites DEL/DUP/INV to symbolic ALT (`<DEL>`/`<DUP>`/`<INV>`).
     Default 50bp. INS is always kept sequence-resolved regardless of this threshold — see design notes.
   - `contigs`: the contig scatter list for stage 2 — set this to match your reference's naming (`chr1` vs `1`, etc).
   - Tune `stage1_bnd_distance` / `stage1_overlap` (intra-sample, cross-caller matching) and `stage2_bnd_distance` / `stage2_overlap` (cross-sample matching) if needed — defaults follow the values discussed earlier (300bp / 0.5 for stage 1, 500bp / 0.5 for stage 2).
   - `max_missing_fraction`: sites with more missing genotypes than this are dropped in the final filter.

4. Dry-run first, then execute:

   ```bash
   snakemake -s Snakefile -n --use-conda               # dry run — check the DAG
   snakemake -s Snakefile -j <ncores> --use-conda       # run locally
   # or, for a cluster: snakemake -s Snakefile --use-conda --profile <your-cluster-profile>
   ```

Output: `results/cohort/cohort.final.vcf.gz`

## Design notes

**Why `svdb --merge` and not Jasmine/Truvari.** All three are reasonable choices for this task. SVDB was used here because its `--data_profile sample` vs `--data_profile cohort` presets map directly onto the two merge tolerances this framework needs (strict same-individual matching vs permissive cross-individual matching), and its `SUPP_VEC`/`FOUNDBY` output gives caller-support provenance for free. If you'd rather use Jasmine or Truvari, swap the `merge_per_sample` and `cohort_merge_contig` rules — the pre/post-processing (normalization, standardization, tagging, tiering, QC) stays the same. Truvari in particular is worth considering if you want sequence-aware benchmarking against a truth set (e.g. GIAB) — see the note on `standardize_sv.py` below, which was designed with truvari compatibility in mind even though svdb is what's wired in here.

**Why sample renaming happens during normalization.** Manta and Dysgu often write a generic or BAM-derived sample name into the VCF's sample column rather than your `sample_id`. Left alone, this causes silent sample collisions during the cohort merge (two different people's calls landing in a column both named e.g. `SAMPLE`). `normalize_caller_vcf` reheaders every VCF to the `sample_id` from `samples.tsv` before anything else happens.

**Why `standardize_sv.py` exists, and why INS is treated differently to DEL/DUP/INV.** Manta and Dysgu represent the same kind of event quite differently in REF/ALT — Dysgu emits minimal, spec-compliant sequence-resolved insertions, while Manta's insertions often carry extra CIGAR-driven padding (e.g. `1M376I4D`), and its large deletions/duplications/inversions are usually already symbolic. Left alone, this skews SVLEN/END comparisons *before* the callers are even merged (`svdb`'s matching relies on SVLEN/END, not REF/ALT sequence), and leaves a visually inconsistent final VCF. `standardize_sv.py` runs after `normalize_caller_vcf` and before the per-sample `svdb` merge, so both callers enter the merge with consistent SVLEN/END and matching gets more accurate, not just cosmetically cleaner output. DEL/DUP/INV at or above `symbolic_threshold` are rewritten to symbolic ALT with recomputed SVLEN/END; INS is always kept sequence-resolved and rewritten to the minimal 1bp-anchor form, regardless of size. That asymmetry is deliberate: truvari (and most SV benchmarking) can reference-resolve symbolic `<DEL>`/`<INV>` (and `<DUP>` via `--dup-to-ins`) for sequence comparison, but it can't recover a novel insertion's sequence from the reference — collapsing INS to symbolic `<INS>` would silently disable sequence-similarity scoring for every insertion.

When Manta supplies a CIGAR, the script parses it to recover the true inserted sequence (CIGAR ops walked over ALT: `M`/`I` consume ALT, `D` does not) rather than trusting the raw REF/ALT as reported. This uncovered a real Manta quirk worth knowing about: some Manta INS records have a CIGAR whose `M`+`I` bases don't actually add up to the length of the reported ALT — the CIGAR and REF/ALT are internally inconsistent in the caller's own output. Rather than guess and risk emitting a silently-truncated insertion sequence, the script detects this mismatch and leaves the record untouched, flagged `INFO/SV_STD_SKIPPED=CIGAR_ALT_MISMATCH`, so those sites are easy to find and inspect. It's worth checking after a real run how common this is in your data.

`standardize_sv.py` also drops BND/TRA records whose *mate* (parsed from BND ALT bracket notation, or `INFO/CHR2` for Dysgu-style records) falls outside the primary contig set (`--primary-contigs`, same list `normalize_caller_vcf` already restricts each record's own CHROM to). `normalize_caller_vcf`'s CHROM-only restriction doesn't catch a BND/TRA whose own CHROM is primary but whose mate points into a decoy/unplaced contig — exactly the kind of alignment artifact decoys exist to absorb. 

**Why `convert_manta_inversions` runs first, and why there's no container involved.** Manta represents inversions as a pair of BND (breakend) records rather than a single `SVTYPE=INV` record. Manta ships its own fix for this, `convertInversion.py`, which needs to run on Manta's *untouched* native output — hence it's Stage -1, before even `normalize_caller_vcf`. Reading the actual script (github.com/Illumina/manta, `src/python/libexec/convertInversion.py`) shows it only imports Python stdlib (`sys`, `gzip`, `io.BufferedReader`, `subprocess.check_output`, `os.path`) and shells out to `samtools` for exactly one thing (`samtools faidx`) — no sibling libexec files, no Manta binaries. That means a local copy of the single file plus a small dedicated conda env (`envs/manta_convert.yml`: python=2.7 + samtools + htslib) covers everything it needs; no Manta container is required at run time, which sidesteps a long list of container-specific problems hit along the way (build-tree files lacking the execute bit inside read-only images, `$MANTA_INSTALL_FOLDER` not being exported, no bundled standalone `samtools` CLI in some Manta images, and Snakemake's `conda:` + `singularity:` combination requiring conda to already be installed *inside* the image, which none of the images tried had).

One real bug in Manta's own script surfaced against actual data: HLA/IMGT decoy contigs in some references use colon-separated allele nomenclature in their names (e.g. `HLA-DRB1*04:03:01`), which breaks
`convertInversion.py`'s assumption that BND ALT bracket notation (`[chrom:pos[`) contains exactly one colon, crashing with `ValueError: too many values to unpack`. Since these are BND calls into decoy MHC-region contigs anyway — the same kind of non-primary-contig-mate artifact `standardize_sv.py` filters out downstream — `convert_manta_inversions` filters them out *before* calling `convertInversion.py` (matching on `*` in ALT, which only appears in these HLA allele names, never in a real contig name), rather than patching Illumina's vendored script.

**Why `FORMAT/CS` instead of relying on SVDB's own `set`/`FOUNDBY` fields.** SVDB re-adds `VARID`/`set`/`FOUNDBY`/`SUPP_VEC` INFO fields at every merge call, and their names collide across stage 1 and stage 2 runs (re-merging a VCF that already has a `FOUNDBY` INFO field produces a duplicate, spec-violating INFO key). `scripts/add_caller_support.py` converts stage 1's `SUPP_VEC` into a per-sample `FORMAT/CS` tag and strips SVDB's provenance INFO fields, so the record is clean going into stage 2. Because `merge_per_sample` uses `--priority 1,2`, a record supported by both callers keeps Manta's REF/ALT/POS/QUAL/FILTER/INFO and genotype/ FORMAT values outright rather than blending the two — Dysgu's specific values for that site aren't retained, only the fact that it also called the event, via `CS=BOTH`. `scripts/tier_cohort.py` does the same provenance cleanup after the cohort merge and additionally aggregates `FORMAT/CS` across all samples into `INFO/CALLER_TIER` (`DUAL_CALLER_MAJORITY` / `MANTA_MAJORITY` / `DYSGU_MAJORITY` / `MIXED` / `NO_CALLS`) plus `N_CALLED`/`N_BOTH`/`N_MANTA_ONLY`/`N_DYSGU_ONLY` counts, so you can filter or stratify downstream analyses by caller support without re-deriving it.

**Why `filter_cohort` checks `FILTER=PASS`, not just `F_MISSING`/`AF`.** svdb's merge carries the FILTER value through from whichever caller's record it kept as the merge representative — that can be a low-confidence call (e.g. Dysgu's `lowProb`, occasionally with `QUAL=0`) that happens to still clear the `F_MISSING`/`AF` thresholds. `filter_cohort` now requires `FILTER="PASS"` (or `.`, since some merge steps leave it unset rather than explicitly `PASS`) alongside the existing thresholds.

**Why contig scattering, not sample scattering, at stage 2.** Stage 1 is already embarrassingly parallel over samples. Stage 2 has to look across all samples at once by construction (that's the point of a cohort merge), so the only axis left to parallelize is genomic — hence scattering by contig via `subset_contig` before `cohort_merge_contig`. 

  What's *not* included: joint/force genotyping

As discussed, after the cohort merge a sample with no record at a site is genuinely missing (`./.`), not confidently reference — Manta or Dysgu simply may not have called anything there, which is different from having evidence of absence. Closing that gap requires re-examining each sample's reads at each merged site, which needs BAM/CRAM access and a genotyper such as Paragraph, GraphTyper2, or SVTyper — out of scope for a VCF-only workflow. If you have BAMs available, the natural place to add this is a new stage between `filter_cohort`'s input and output: take `cohort.annotated.vcf.gz`, force-genotype every `./.` cell against the corresponding sample's BAM, then re-run the AF/F_MISSING annotation and filter on the result. 
## What was validated

The full DAG (Manta inversion conversion → normalization → SV standardization → per-sample merge → tagging → contig scatter → cohort merge → tiering → concat → annotate → filter) was run end-to-end against small synthetic Manta/Dysgu VCFs, including a multi-sample case, and produced a correctly-tiered multi-sample VCF with accurate `CALLER_TIER`/`N_BOTH`/`N_MANTA_ONLY`/`N_DYSGU_ONLY`/`AF`/`F_MISSING` annotations. 

`standardize_sv.py`'s core logic (CIGAR parsing, insertion-sequence extraction, CIGAR/ALT consistency checking, mate-contig extraction from both BND bracket notation and `INFO/CHR2`, END-filling for symbolic INS/TRA without a native END) was unit-tested directly against real Manta and Dysgu example records rather than synthetic data — including the CIGAR/ALT-mismatch case above, found this way. 

`convert_manta_inversions` was debugged iteratively against real HG002 Manta output on the target cluster; issues found and fixed along the way included pysam's reserved `INFO/END` field needing `.stop` instead of `info["END"]`, pysam returning tuples for `Number=.`/`Number=A` INFO fields like Manta's `CIGAR` and `CHR2`, several container-specific permission/PATH/samtools-discovery problems (resolved by dropping the container in favor of a local script copy + dedicated conda env — see design notes), a `set -u`/`set -e` bash-strict-mode trap from Snakemake's shell templating, and the HLA-decoy contig-name bug described above. As of the last run, it was working through real data cleanly with the HLA filter in place; a full end-to-end run against the real 234-sample cohort has not yet been confirmed. Before a production run, worth checking: 

- Manta emits `BND`-notation records for some deletions/duplications in addition to `<DEL>`/`<DUP>` symbolic ALTs; confirm SVDB's BND matching (`--bnd_distance`) behaves as expected on your actual VCFs.
- At 234 samples, the `--vcf` argument list passed to `cohort_merge_contig` will be long (234 file paths per contig); this worked in testing but if you hit shell argument-length limits, switch that rule to write the file list to a temp file (SVDB does not   currently support `--vcf @filelist`, so file-list support would need to be added or file paths shortened/symlinked into a shorter common directory).
- Tune `stage1_bnd_distance`/`overlap` and `stage2_bnd_distance`/`overlap` against a small truth set if you have one (e.g. a few samples with orthogonal validation) rather than trusting the defaults blindly.
- Check how many records `standardize_sv.py` flags with `SV_STD_SKIPPED=CIGAR_ALT_MISMATCH` and how many `convert_manta_inversions` drops for HLA-decoy ALTs (both print counts to stderr) — high counts of either are worth a closer look at your specific Manta build/reference combination.