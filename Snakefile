"""
Manta + Dysgu structural-variant merge workflow (cohort scale, VCF-only).
 
Stage -1 (Manta only): run Manta's own convertInversion.py (a local copy,
    Python 2 via a small dedicated conda env) on the raw Manta VCF first,
    converting its BND-pair representation of inversions into proper
    SVTYPE=INV records. Must run before normalize_caller_vcf, on Manta's
    untouched native output.
 
Stage 1 (per sample, parallel over all samples):
    normalize Manta + Dysgu VCFs -> standardize SV representation
    -> svdb merge (data_profile=sample) -> tag FORMAT/CS caller support
    -> sort/index
 
Stage 2 (per contig, parallel over contigs):
    subset every sample's stage-1 VCF to the contig -> svdb merge across
    samples (data_profile=cohort) -> tier caller support at cohort level
    -> sort/index
 
Stage 3:
    concat contigs
    -> [if config["enable_gcnv"]: reconcile DEL/DUP against gCNV depth
       calls via GATK SVCluster -- see "Stage 3.5" below]
    -> annotate AF / F_MISSING -> filter -> final cohort VCF
 
Stage 3.5 (optional, config["enable_gcnv"], depth-caller integration):
    gCNV calling itself runs outside this pipeline (cohort-mode
    GermlineCNVCaller needs BAM/CRAM access and is jointly across all
    samples at once, not per-sample like Manta/Dysgu). Given each
    sample's raw PostprocessGermlineCNVCalls VCF: build a ploidy table
    from samples.tsv's "sex" column -> standardize each gCNV VCF (add
    INFO/ALGORITHMS, FORMAT/ECN) -> split the cohort VCF into DEL/DUP vs
    everything else (gCNV can only corroborate copy-number types) ->
    GATK SVCluster reconciles the DEL/DUP split against all samples'
    standardized gCNV VCFs (svdb has no per-pair tolerance support, so it
    can't distinguish depth-caller breakpoint fuzziness from Manta/
    Dysgu's tighter PESR precision the way SVCluster natively does) ->
    recombine with the untouched non-DEL/DUP records.
 
Joint/force-genotyping (filling in samples with no caller support at a
site with an evidence-based genotype call) is NOT included here: it
requires read-level access to each sample's BAM/CRAM via a tool such as
Paragraph or GraphTyper2, which is out of scope for a VCF-only workflow.
See README.md for notes on adding that step.
"""
 
import pandas as pd
 
configfile: "/software/projects/pawsey0933/SVmerge/config.yaml"
 
SAMPLES_DF = pd.read_csv(config["samples_tsv"], sep="\t", dtype=str).set_index(
    "sample_id", drop=False
)
SAMPLES = SAMPLES_DF["sample_id"].tolist()
CONTIGS = config["contigs"]
OUTDIR = config["outdir"]
 
 
def _manta_vcf(sample):
    return SAMPLES_DF.loc[sample, "manta_vcf"]
 
 
def _dysgu_vcf(sample):
    return SAMPLES_DF.loc[sample, "dysgu_vcf"]
 
 
def _gcnv_vcf(sample):
    return SAMPLES_DF.loc[sample, "gcnv_vcf"]
 
 
# gCNV/SVCluster depth-caller integration is opt-in (config["enable_gcnv"]):
# it needs samples.tsv to have "gcnv_vcf" and "sex" columns that a plain
# Manta+Dysgu run doesn't require, so the extra rules below are only
# pulled into the DAG -- and _gcnv_vcf()/the "sex" column only looked up
# -- when a rule that's actually part of the requested output references
# them. See annotate_af_missingness's input for the branch point.
ENABLE_GCNV = config.get("enable_gcnv", False)
 
# Ploidy table source (see config.yaml's "gcnv_ploidy_calls_dir" comment):
# by default we derive it from samples.tsv's "sex" column
# (rule make_ploidy_table); if gcnv_ploidy_calls_dir is set, we instead
# reuse DetermineGermlineContigPloidy's own per-sample calls
# (rule ploidy_table_from_gcnv_calls). Only one of the two rules actually
# runs for a given config, selected via _ploidy_table() below -- both are
# always present in the DAG description, but Snakemake only builds
# whichever one the config's chosen output path resolves to.
GCNV_PLOIDY_CALLS_DIR = config.get("gcnv_ploidy_calls_dir")
 
 
def _ploidy_table(wc):
    if GCNV_PLOIDY_CALLS_DIR:
        return f"{OUTDIR}/svcluster/ploidy_table.from_gcnv_calls.tsv"
    return f"{OUTDIR}/svcluster/ploidy_table.from_sex_column.tsv"
 
 
rule all:
    input:
        f"{OUTDIR}/cohort/cohort.final.vcf.gz",
        f"{OUTDIR}/cohort/cohort.final.vcf.gz.tbi",
 
 
# --------------------------------------------------------------------------
# Stage -1 (Manta only): convert Manta's BND-pair inversion representation
# to SVTYPE=INV via Manta's own convertInversion.py. Run on the raw Manta
# VCF, before any of the generic sort/reheader/norm steps below.
#
# No container here: convertInversion.py (read directly from Illumina's
# repo -- github.com/Illumina/manta/blob/master/src/python/libexec/
# convertInversion.py) only imports Python stdlib (sys, gzip,
# io.BufferedReader, subprocess.check_output, os.path) and shells out to
# samtools for exactly one thing (`samtools faidx`), so a local copy of
# just that one file plus a small dedicated conda env covers everything
# it needs -- no sibling libexec files, no Manta binaries, no container.
# It genuinely does need Python 2 though (not just "written for" it): it
# reads the VCF in binary mode and does str/bytes operations on the lines
# that raise TypeError under Python 3.
#
# One-time setup (see envs/manta_convert.yml for the env):
#   singularity exec /software/projects/pawsey0933/manta/manta_main.sif \
#       cat /build/manta/build/src/python/libexec/convertInversion.py \
#       > <config["manta_convertinversion_script"]>
# --------------------------------------------------------------------------
 
rule convert_manta_inversions:
    input:
        vcf=lambda wc: _manta_vcf(wc.sample),
    output:
        vcf=f"{OUTDIR}/manta_inv_converted/{{sample}}.manta.invconv.vcf.gz",
        tbi=f"{OUTDIR}/manta_inv_converted/{{sample}}.manta.invconv.vcf.gz.tbi",
    params:
        ref=config["reference"],
        tmp=lambda wc: f"{OUTDIR}/manta_inv_converted/{wc.sample}.manta.invconv.vcf",
    conda: "envs/manta_convert.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/manta_inv_converted

        # envs/manta_convert.yml's python is python2.7, so plain "python"
        # unambiguously means the right interpreter here.
        SAMTOOLS_BIN=$(command -v samtools || true)
        if [ -z "$SAMTOOLS_BIN" ]; then
            echo "ERROR: samtools not found on PATH -- envs/manta_convert.yml" \
                 "should have installed it into this conda env." >&2
            exit 1
        fi

        # convertInversion.py assumes BND ALT bracket notation
        # ([chrom:pos[ or ]chrom:pos]) has exactly one colon. HLA/IMGT
        # decoy contigs in this reference use colon-separated allele
        # nomenclature themselves (e.g. HLA-DRB1*04:03:01), which breaks
        # that assumption and crashes the script with "ValueError: too
        # many values to unpack". These are BND calls into decoy
        # MHC-region contigs -- exactly the kind of non-primary-contig
        # mate that standardize_sv.py already drops downstream -- so
        # filter them out here before convertInversion.py ever sees them,
        # rather than patch Illumina's vendored script. "*" only appears
        # in these HLA allele names, never in a real contig name, so it's
        # a reliable marker.

        FILTERED_VCF="{params.tmp}.hla_filtered.vcf.gz"
        zcat {input.vcf} | awk -F'\t' '$1 ~ /^#/ || $5 !~ /\*/' | bgzip > "$FILTERED_VCF"
        N_TOTAL=$(zcat {input.vcf} | grep -vc '^#' || true)
        N_KEPT=$(zcat "$FILTERED_VCF" | grep -vc '^#' || true)
        echo "[convert_manta_inversions] dropped $((N_TOTAL - N_KEPT)) / $N_TOTAL records" \
             "with an HLA-decoy-contig ALT (unparseable by convertInversion.py)" >&2

        python scripts/convertInversion.py "$SAMTOOLS_BIN" {params.ref} "$FILTERED_VCF" \
            > {params.tmp}
        bgzip -f {params.tmp}
        tabix -f -p vcf {output.vcf}
        """


# --------------------------------------------------------------------------
# Stage 0: normalize each sample's caller VCFs
# --------------------------------------------------------------------------
 
rule normalize_caller_vcf:
    input:
        vcf=lambda wc: (
            f"{OUTDIR}/manta_inv_converted/{wc.sample}.manta.invconv.vcf.gz"
            if wc.caller == "manta"
            else _dysgu_vcf(wc.sample)
        ),
    output:
        vcf=f"{OUTDIR}/normalized/{{sample}}.{{caller}}.norm.vcf.gz",
        tbi=f"{OUTDIR}/normalized/{{sample}}.{{caller}}.norm.vcf.gz.tbi",
    params:
        ref=config["reference"],
        tmp=lambda wc: f"{OUTDIR}/normalized/{wc.sample}.{wc.caller}.sorted.tmp.vcf.gz",
    wildcard_constraints:
        caller="manta|dysgu",
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/normalized
        bcftools sort -Oz -o {params.tmp} {input.vcf}
        # Rename the (single) sample column to this sample's sample_id.
        # Manta/Dysgu often write a generic or BAM-derived sample name;
        # without this, distinct samples can collide onto the same column
        # name once merged, silently dropping samples.
        printf '%s\n' {wildcards.sample} > {params.tmp}.rename.txt
        bcftools reheader -s {params.tmp}.rename.txt {params.tmp} -o {params.tmp}.renamed.vcf.gz
        # Restrict to primary assembly contigs (chr1-22,X,Y). Dysgu emits spurious
        # calls at/past the end of minor contigs (decoys, _random/_alt scaffolds,
        # chrM) -- e.g. chrM:16570 (chrM is only 16569bp) and
        # chr9_KI270718v1_random:38055 -- which trip bcftools norm's reference
        # check. These contigs aren't in the GIAB benchmark regions anyway, so
        # derive the primary contig list from the reference's .fai rather than
        # excluding problem contigs one at a time.
        PRIMARY_CONTIGS=$(awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {{print $1}}' {params.ref}.fai | paste -sd, -)
        bcftools view -t "${{PRIMARY_CONTIGS}}" -Oz -o {params.tmp}.primary.vcf.gz {params.tmp}.renamed.vcf.gz
        bcftools norm -f {params.ref} -m -any -Oz -o {output.vcf} {params.tmp}.primary.vcf.gz
        tabix -f -p vcf {output.vcf}
        rm -f {params.tmp} {params.tmp}.rename.txt {params.tmp}.renamed.vcf.gz {params.tmp}.primary.vcf.gz
        """
 
 
# --------------------------------------------------------------------------
# Stage 0.5: harmonize SV representation across callers (see
# scripts/standardize_sv.py docstring for rationale). Runs before the
# per-sample svdb merge so that Manta/Dysgu records compare like-for-like
# SVLEN/END, not just so the final VCF looks uniform.
# --------------------------------------------------------------------------
 
rule standardize_sv:
    input:
        vcf=f"{OUTDIR}/normalized/{{sample}}.{{caller}}.norm.vcf.gz",
        tbi=f"{OUTDIR}/normalized/{{sample}}.{{caller}}.norm.vcf.gz.tbi",
    output:
        vcf=f"{OUTDIR}/standardized/{{sample}}.{{caller}}.std.vcf.gz",
        tbi=f"{OUTDIR}/standardized/{{sample}}.{{caller}}.std.vcf.gz.tbi",
    params:
        tmp=lambda wc: f"{OUTDIR}/standardized/{wc.sample}.{wc.caller}.std.tmp.vcf",
        symbolic_threshold=config.get("symbolic_threshold", 50),
        ref=config["reference"],
    wildcard_constraints:
        caller="manta|dysgu",
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/standardized
        # Same primary-contig set normalize_caller_vcf already restricted
        # record CHROM to. Passed through here so standardize_sv.py can
        # also drop BND/TRA records whose *mate* (not just their own
        # CHROM) points at a decoy/unplaced contig -- CHROM-only filtering
        # doesn't catch those, since the record's own CHROM is primary.
        PRIMARY_CONTIGS=$(awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {{print $1}}' {params.ref}.fai | paste -sd, -)
        python3 scripts/standardize_sv.py {input.vcf} {params.tmp} \
            --symbolic-threshold {params.symbolic_threshold} \
            --primary-contigs "${{PRIMARY_CONTIGS}}"
        bcftools sort -Oz -o {output.vcf} {params.tmp}
        tabix -f -p vcf {output.vcf}
        rm -f {params.tmp}
        """
 
 
# --------------------------------------------------------------------------
# Stage 1: per-sample Manta + Dysgu merge
# --------------------------------------------------------------------------
 
rule merge_per_sample:
    input:
        manta=f"{OUTDIR}/standardized/{{sample}}.manta.std.vcf.gz",
        dysgu=f"{OUTDIR}/standardized/{{sample}}.dysgu.std.vcf.gz",
    output:
        vcf=temp(f"{OUTDIR}/per_sample_merged/{{sample}}.raw_merge.vcf"),
    params:
        max_dist=config["stage1_bnd_distance"],
        overlap=config["stage1_overlap"],
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/per_sample_merged
        svdb --merge \
            --vcf {input.manta}:1 {input.dysgu}:2 \
            --priority 1,2 \
            --bnd_distance {params.max_dist} \
            --overlap {params.overlap} \
            --data_profile sample \
            > {output.vcf}
        """
 
 
rule tag_caller_support:
    input:
        f"{OUTDIR}/per_sample_merged/{{sample}}.raw_merge.vcf",
    output:
        vcf=temp(f"{OUTDIR}/per_sample_merged/{{sample}}.tagged.vcf"),
    conda: "envs/svmerge.yml"
    shell:
        r"""
        python3 scripts/add_caller_support.py {input} {output.vcf}
        """
 
 
rule sort_index_per_sample:
    input:
        f"{OUTDIR}/per_sample_merged/{{sample}}.tagged.vcf",
    output:
        vcf=f"{OUTDIR}/per_sample_merged/{{sample}}.merged.vcf.gz",
        tbi=f"{OUTDIR}/per_sample_merged/{{sample}}.merged.vcf.gz.tbi",
    conda: "envs/svmerge.yml"
    shell:
        r"""
        bcftools sort -Oz -o {output.vcf} {input}
        tabix -f -p vcf {output.vcf}
        """
 
 
# --------------------------------------------------------------------------
# Stage 2: cohort-level cross-sample merge, scattered by contig
# --------------------------------------------------------------------------
 
rule subset_contig:
    input:
        vcf=f"{OUTDIR}/per_sample_merged/{{sample}}.merged.vcf.gz",
        tbi=f"{OUTDIR}/per_sample_merged/{{sample}}.merged.vcf.gz.tbi",
    output:
        vcf=temp(f"{OUTDIR}/per_sample_by_contig/{{sample}}.{{contig}}.vcf.gz"),
        tbi=temp(f"{OUTDIR}/per_sample_by_contig/{{sample}}.{{contig}}.vcf.gz.tbi"),
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/per_sample_by_contig
        bcftools view -r {wildcards.contig} -Oz -o {output.vcf} {input.vcf}
        tabix -f -p vcf {output.vcf}
        """
 
 
rule cohort_merge_contig:
    input:
        vcfs=lambda wc: expand(
            f"{OUTDIR}/per_sample_by_contig/{{sample}}.{wc.contig}.vcf.gz",
            sample=SAMPLES,
        ),
    output:
        vcf=temp(f"{OUTDIR}/cohort_by_contig/{{contig}}.raw_merge.vcf"),
    params:
        max_dist=config["stage2_bnd_distance"],
        overlap=config["stage2_overlap"],
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/cohort_by_contig
        svdb --merge --no_tag \
            --vcf {input.vcfs} \
            --bnd_distance {params.max_dist} \
            --overlap {params.overlap} \
            --data_profile cohort \
            > {output.vcf}
        """
 
 
rule tier_contig:
    input:
        f"{OUTDIR}/cohort_by_contig/{{contig}}.raw_merge.vcf",
    output:
        vcf=temp(f"{OUTDIR}/cohort_by_contig/{{contig}}.tiered.vcf"),
    conda: "envs/svmerge.yml"
    shell:
        r"""
        python3 scripts/tier_cohort.py {input} {output.vcf}
        """
 
 
rule sort_index_contig:
    input:
        f"{OUTDIR}/cohort_by_contig/{{contig}}.tiered.vcf",
    output:
        vcf=f"{OUTDIR}/cohort_by_contig/{{contig}}.tiered.vcf.gz",
        tbi=f"{OUTDIR}/cohort_by_contig/{{contig}}.tiered.vcf.gz.tbi",
    conda: "envs/svmerge.yml"
    shell:
        r"""
        bcftools sort -Oz -o {output.vcf} {input}
        tabix -f -p vcf {output.vcf}
        """
 
 
# --------------------------------------------------------------------------
# Stage 3.5 (optional, config["enable_gcnv"]): prep for depth-caller
# (gCNV) integration -- ploidy table + per-sample gCNV standardization.
# Independent of stages 1-2, so these can build any time; placed here
# since nothing else in the pipeline uses them.
# --------------------------------------------------------------------------
 
rule make_ploidy_table:
    input:
        samples_tsv=config["samples_tsv"],
    output:
        tsv=f"{OUTDIR}/svcluster/ploidy_table.from_sex_column.tsv",
    params:
        contigs=",".join(CONTIGS),
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/svcluster
        python3 scripts/make_ploidy_table.py {input.samples_tsv} "{params.contigs}" {output.tsv}
        """
 
 
# Alternative to make_ploidy_table above: reuses DetermineGermlineContigPloidy's
# own per-sample ploidy calls (from Snakefile.gcnv's determine_germline_contig_ploidy
# rule) instead of deriving ploidy from samples.tsv's "sex" column. See
# config.yaml's "gcnv_ploidy_calls_dir" comment and
# scripts/ploidy_calls_to_table.py's docstring for why this is preferred
# when available. Only pulled into the DAG when gcnv_ploidy_calls_dir is set
# -- see _ploidy_table() above for the branch point.
rule ploidy_table_from_gcnv_calls:
    input:
        ploidy_calls_dir=directory(GCNV_PLOIDY_CALLS_DIR) if GCNV_PLOIDY_CALLS_DIR else [],
    output:
        tsv=f"{OUTDIR}/svcluster/ploidy_table.from_gcnv_calls.tsv",
    params:
        contigs=",".join(CONTIGS),
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/svcluster
        python3 scripts/ploidy_calls_to_table.py {input.ploidy_calls_dir} "{params.contigs}" {output.tsv}
        """
 
 
rule standardize_gcnv:
    input:
        vcf=lambda wc: _gcnv_vcf(wc.sample),
        ploidy=_ploidy_table,
    output:
        vcf=f"{OUTDIR}/gcnv_standardized/{{sample}}.gcnv.std.vcf.gz",
        tbi=f"{OUTDIR}/gcnv_standardized/{{sample}}.gcnv.std.vcf.gz.tbi",
    params:
        tmp=lambda wc: f"{OUTDIR}/gcnv_standardized/{wc.sample}.gcnv.std.tmp.vcf",
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/gcnv_standardized
        python3 scripts/standardize_gcnv.py {input.vcf} {params.tmp} {input.ploidy} {wildcards.sample}
        bcftools sort -Oz -o {output.vcf} {params.tmp}
        tabix -f -p vcf {output.vcf}
        rm -f {params.tmp}
        """
 
 
# --------------------------------------------------------------------------
# Stage 3: concat contigs, annotate, QC filter
# --------------------------------------------------------------------------
 
rule concat_contigs:
    input:
        vcfs=expand(f"{OUTDIR}/cohort_by_contig/{{contig}}.tiered.vcf.gz", contig=CONTIGS),
        tbis=expand(f"{OUTDIR}/cohort_by_contig/{{contig}}.tiered.vcf.gz.tbi", contig=CONTIGS),
    output:
        vcf=f"{OUTDIR}/cohort/cohort.merged.vcf.gz",
        tbi=f"{OUTDIR}/cohort/cohort.merged.vcf.gz.tbi",
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/cohort
        bcftools concat -Oz -o {output.vcf} {input.vcfs}
        tabix -f -p vcf {output.vcf}
        """
 
 
# --------------------------------------------------------------------------
# Stage 3.5 continued (optional, config["enable_gcnv"]): reconcile the
# concatenated cohort VCF's DEL/DUP records against the standardized gCNV
# VCFs via GATK SVCluster, then recombine with the untouched non-DEL/DUP
# records. Only pulled into the DAG when config["enable_gcnv"] is set --
# see annotate_af_missingness's input below for the branch point.
# --------------------------------------------------------------------------
 
rule split_deldup:
    input:
        vcf=f"{OUTDIR}/cohort/cohort.merged.vcf.gz",
        tbi=f"{OUTDIR}/cohort/cohort.merged.vcf.gz.tbi",
    output:
        deldup=f"{OUTDIR}/svcluster/cohort.deldup.vcf.gz",
        deldup_tbi=f"{OUTDIR}/svcluster/cohort.deldup.vcf.gz.tbi",
        rest=f"{OUTDIR}/svcluster/cohort.rest.vcf.gz",
        rest_tbi=f"{OUTDIR}/svcluster/cohort.rest.vcf.gz.tbi",
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/svcluster
        # gCNV is a depth-based caller: it can only corroborate copy-number
        # types (DEL/DUP), so only those go through SVCluster reconciliation.
        # Everything else (INS/INV/BND/...) passes through untouched.
        bcftools view -i 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP"' \
            -Oz -o {output.deldup} {input.vcf}
        tabix -f -p vcf {output.deldup}
        bcftools view -e 'INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP"' \
            -Oz -o {output.rest} {input.vcf}
        tabix -f -p vcf {output.rest}
        """
 
 
rule svcluster_depth_merge:
    input:
        deldup=f"{OUTDIR}/svcluster/cohort.deldup.vcf.gz",
        deldup_tbi=f"{OUTDIR}/svcluster/cohort.deldup.vcf.gz.tbi",
        gcnv_vcfs=expand(
            f"{OUTDIR}/gcnv_standardized/{{sample}}.gcnv.std.vcf.gz", sample=SAMPLES
        ),
        gcnv_tbis=expand(
            f"{OUTDIR}/gcnv_standardized/{{sample}}.gcnv.std.vcf.gz.tbi", sample=SAMPLES
        ),
        ploidy=_ploidy_table,
    output:
        vcf=f"{OUTDIR}/svcluster/deldup.clustered.vcf.gz",
        tbi=f"{OUTDIR}/svcluster/deldup.clustered.vcf.gz.tbi",
    params:
        ref=config["reference"],
        algorithm=config.get("svcluster_algorithm", "SINGLE_LINKAGE"),
        gcnv_args=lambda wc, input: " ".join(f"-V {v}" for v in input.gcnv_vcfs),
    conda: "envs/gatk_svcluster.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/svcluster
        REF="{params.ref}"
        REF_DICT="${{REF%.*}}.dict"
        if [ ! -f "$REF_DICT" ]; then
            gatk CreateSequenceDictionary -R "$REF" -O "$REF_DICT"
        fi
        gatk SVCluster \
            -V {input.deldup} \
            {params.gcnv_args} \
            --ploidy-table {input.ploidy} \
            -R "$REF" \
            --algorithm {params.algorithm} \
            -O {output.vcf}
        tabix -f -p vcf {output.vcf}
        """
 
 
rule recombine_deldup:
    input:
        clustered=f"{OUTDIR}/svcluster/deldup.clustered.vcf.gz",
        clustered_tbi=f"{OUTDIR}/svcluster/deldup.clustered.vcf.gz.tbi",
        rest=f"{OUTDIR}/svcluster/cohort.rest.vcf.gz",
        rest_tbi=f"{OUTDIR}/svcluster/cohort.rest.vcf.gz.tbi",
    output:
        vcf=f"{OUTDIR}/cohort/cohort.merged.with_depth.vcf.gz",
        tbi=f"{OUTDIR}/cohort/cohort.merged.with_depth.vcf.gz.tbi",
    conda: "envs/svmerge.yml"
    shell:
        r"""
        mkdir -p {OUTDIR}/cohort
        bcftools concat -a {input.clustered} {input.rest} -Oz -o {output.vcf}.tmp.gz
        bcftools sort -Oz -o {output.vcf} {output.vcf}.tmp.gz
        tabix -f -p vcf {output.vcf}
        rm -f {output.vcf}.tmp.gz
        """
 
 
def _cohort_merged_vcf(wc):
    if ENABLE_GCNV:
        return f"{OUTDIR}/cohort/cohort.merged.with_depth.vcf.gz"
    return f"{OUTDIR}/cohort/cohort.merged.vcf.gz"
 
 
def _cohort_merged_tbi(wc):
    return _cohort_merged_vcf(wc) + ".tbi"
 
 
rule annotate_af_missingness:
    input:
        vcf=_cohort_merged_vcf,
        tbi=_cohort_merged_tbi,
    output:
        vcf=f"{OUTDIR}/cohort/cohort.annotated.vcf.gz",
        tbi=f"{OUTDIR}/cohort/cohort.annotated.vcf.gz.tbi",
    conda: "envs/svmerge.yml"
    shell:
        r"""
        bcftools +fill-tags {input.vcf} -Oz -o {output.vcf} -- -t AF,F_MISSING
        tabix -f -p vcf {output.vcf}
        """
 
 
rule filter_cohort:
    input:
        vcf=f"{OUTDIR}/cohort/cohort.annotated.vcf.gz",
        tbi=f"{OUTDIR}/cohort/cohort.annotated.vcf.gz.tbi",
    output:
        vcf=f"{OUTDIR}/cohort/cohort.final.vcf.gz",
        tbi=f"{OUTDIR}/cohort/cohort.final.vcf.gz.tbi",
    params:
        max_missing=config["max_missing_fraction"],
        min_maf=config["min_maf"],
    conda: "envs/svmerge.yml"
    shell:
        r"""
        # FILTER="PASS" excludes low-confidence calls (e.g. Dysgu's
        # "lowProb") that svdb otherwise carries through untouched from
        # whichever input record it kept as the merge representative.
        # FILTER="." is also accepted since some callers/merge steps
        # leave FILTER unset rather than explicitly PASS.
        bcftools view \
            -i '(FILTER="PASS" || FILTER=".") && INFO/F_MISSING<={params.max_missing} && (INFO/AF>={params.min_maf} && INFO/AF<=(1-{params.min_maf}))' \
            -Oz -o {output.vcf} {input.vcf}
        tabix -f -p vcf {output.vcf}
        """