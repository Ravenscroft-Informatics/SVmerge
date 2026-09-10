# SVmerge

A Snakemake workflow implementing the two-stage merge framework discussed previously: per-sample Manta+Dysgu merge, then cohort-level cross-sample merge, with caller-support tiering carried through both stages. Manta's raw output is first passed through an inversion-representation fix and an SV-representation standardization pass before any merging happens. An optional Stage 3.5 can additionally reconcile the cohort's DEL/DUP calls against depth-based GATK gCNV calls via GATK SVCluster. 

## What this does

```
Manta only:
  convertInversion.py (local copy of Manta's own script) converts Manta's BND-pair representation of inversions into SVTYPE=INV records

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
  concat contigs
    -> [if enable_gcnv: reconcile DEL/DUP against gCNV depth calls via GATK SVCluster, recombine with the untouched non-DEL/DUP records -- see "Stage 3.5" below]
    -> annotate AF / F_MISSING
    -> filter (FILTER=PASS, F_MISSING, AF) -> cohort.final.vcf.gz
```

**Stage 3.5 (optional, enable_gcnv: true): depth-caller integration** 
- gCNV calling itself runs OUTSIDE this pipeline since it needs BAM/CRAM access and calls jointly across the whole cohort at once, not per-sample like `Manta`/`Dysgu`.
- Use (`cnv_snakemake`)[https://github.com/RAVING-Informatics/cnv_snakemake] which performs`CollectReadCounts` -> `DetermineGermlineContigPloidy` -> `GermlineCNVCaller` -> Postprocess `GermlineCNVCalls`

Given each sample's raw `PostprocessGermlineCNVCalls` VCF:
  build a ploidy table from samples.tsv's "sex" column
    -> standardize each sample's gCNV VCF (add `INFO/ALGORITHMS=depth`, `FORMAT/ECN`)
    -> split the concatenated cohort VCF into `DEL`/`DUP` vs everything else (gCNV can only corroborate copy-number types)
    -> `GATK SVCluster` reconciles the `DEL`/`DUP` split against all samples' standardized gCNV VCFs
    -> recombine with the untouched non-`DEL`/`DUP` records

## Requirements

Install `snakemake` using `conda`:

```bash
conda create -n svmerge -c bioconda snakemake=8.30.0
```
Three conda envs are used, activated per-rule via Snakemake's `conda:` directive (run with `snakemake --use-conda ...`):

- `envs/svmerge.yml` — Python 3 + pysam + bcftools/tabix. Used by essentially every rule (`standardize_sv.py`, `add_caller_support.py`, `tier_cohort.py`, all the `bcftools`/`svdb` steps).
- `envs/manta_convert.yml` — Python 2.7 + samtools + htslib, used *only* by `convert_manta_inversions`. Kept separate because Manta's own `convertInversion.py` is Python 2 only (confirmed by reading the script: it does byte/str operations that raise `TypeError` under  Python 3), and conda envs pin a single Python major version.
- `envs/gatk_svcluster.yml` - GATK4, used only by `svcluster_depth_merge` `when enable_gcnv: true`. A JVM tool with no relation to either of the above, so it gets its own minimal env rather than being folded into `svmerge.yml`.

## Setup

1. Fill in `samples.tsv` (see `samples.tsv.example`):

   ```
   sample_id	manta_vcf	dysgu_vcf
   sample001	/path/to/sample001.manta.vcf.gz	/path/to/sample001.dysgu.vcf.gz
   ...
   ```

If you're using Stage 3.5 (`enable_gcnv: true`), two more columns are required:
  ```
    sample_id	manta_vcf	dysgu_vcf	gcnv_vcf	sex
    sample001	.../sample001.manta.vcf.gz	.../sample001.dysgu.vcf.gz	.../sample001.postprocess.vcf.gz	F
    ...
  ```
- `gcnv_vcf`: each sample's raw `PostprocessGermlineCNVCalls` output (from a gCNV run you've done separately — see "Setting up gCNV" below). Not GATK-SV-standardized; `standardize_gcnv.py` handles that here.
- `sex`: `M`/`F` (or `male`/`female`, case-insensitive). Used by `make_ploidy_table.py` to build the `--ploidy-table SVCluster` needs (autosomes always ploidy 2; X/Y derived from this column). This is a fallback, not the preferred source — see `gcnv_ploidy_calls_dir` below.

2. Edit `config.yaml`:
   - `reference`: path to the reference FASTA used for both callers (needs a `.fai` alongside it — `samtools faidx reference.fasta`).
   - `symbolic_threshold`: size (bp) at or above which `standardize_sv.py` rewrites DEL/DUP/INV to symbolic ALT (`<DEL>`/`<DUP>`/`<INV>`). Default 50bp. INS is always kept sequence-resolved regardless of this threshold — see design notes.
   - `contigs`: the contig scatter list for stage 2 — set this to match your reference's naming (`chr1` vs `1`, etc).
   - Tune `stage1_bnd_distance` / `stage1_overlap` (intra-sample, cross-caller matching) and `stage2_bnd_distance` / `stage2_overlap` (cross-sample matching) if needed — defaults follow the values discussed earlier (300bp / 0.5 for stage 1, 500bp / 0.5 for stage 2).
   - `max_missing_fraction`: sites with more missing genotypes than this are dropped in the final filter.
   - `enable_gcnv`: set `true` to turn on Stage 3.5 depth-caller integration (requires `gcnv_vcf`/`sex` columns in `samples.tsv`, per above). Defaults to `false` — a plain Manta+Dysgu run needs no changes here.
   - `svcluster_algorithm`: GATK SVCluster's clustering algorithm (`SINGLE_LINKAGE` default, or `MAX_CLIQUE`/`DEFRAGMENT_CNV` — see the comment above this key in `config.yaml`).
   - `gcnv_ploidy_calls_dir`: optional. If you ran `cnv_snakemake`, `DetermineGermlineContigPloidy` already computed each sample's actual per-contig ploidy from its own coverage — and `GermlineCNVCaller`/ `PostprocessGermlineCNVCalls` used exactly those calls as the baseline copy number when generating the gCNV VCFs. Point this at that run's ploidy-calls directory (`config_gcnv.yaml`'s outdir + `/ploidy/ploidy-calls`) to build the SVCluster ploidy table from those real calls instead of guessing from `samples.tsv`'s sex column — via `scripts/ploidy_calls_to_table.py` rather than `make_ploidy_table.py`. Preferred whenever you have it: it can't miss a real sex-chromosome aneuploidy or a wrong sex value the way the sex-column derivation can, and it guarantees the ploidy table agrees with what the gCNV calls were actually made under. Leave unset (`null`) to keep using the sex-column table.

3. (optional) Run `GATK-gCNV` workflow using (`cnv_snakemake`)[https://github.com/RAVING-Informatics/cnv_snakemake]
4. Dry-run first, then execute:

   ```bash
   snakemake -s Snakefile -n --use-conda               # dry run — check the DAG
   snakemake -s Snakefile -j <ncores> --use-conda       # run locally
   # or, for a cluster: snakemake -s Snakefile --use-conda --profile <your-cluster-profile>
   ```

## Outputs
- `results/cohort/cohort.final.vcf.gz` 

### Filter options
Use bcftools to filter for variants with support from multiple callers:

```
# calls with support from both PESR callers
$ bcftools view -H -i 'INFO/N_BOTH==1' cohort.final.vcf.gz > cohort.final.both.vcf.gz

# calls with support from both PESR callers and depth 
$ bcftools view -H -i 'INFO/ALGORITHMS=="depth" && INFO/ALGORITHMS=="manta" && INFO/ALGORITHMS=="dysgu"' cohort.final.vcf.gz > cohort.final.both.depth.vcf.gz

# calls with support from at least one PESR and depth
$ bcftools view -H -i 'INFO/ALGORITHMS=="depth" && (INFO/ALGORITHMS=="manta" || INFO/ALGORITHMS=="dysgu")' cohort.final.pesr.depth.vcf.gz 
```

## Running on cluster
- Use cluster profile config: `cluster_config.yaml`