#!/bin/bash -l

#SBATCH --job-name=sv_merge
#SBATCH --account=pawsey0933
#SBATCH --partition=debug
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --mail-user=chiara.folland@perkins.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out
#SBATCH --export=ALL

module load singularity/4.1.0-slurm 

#Load conda env
source /software/projects/pawsey0933/cfolland/miniconda3/etc/profile.d/conda.sh
mamba activate /software/projects/pawsey0933/cfolland/miniforge3/envs/drop_env

#Change cache dir
export XDG_CACHE_HOME=/scratch/pawsey0933/cfolland/rna_vc/.cache
export SINGULARITY_CACHEDIR=/scratch/pawsey0933/cfolland/singularity_cache

#Change to working dir
cd /software/projects/pawsey0933/SVmerge

#Unlock working directory
snakemake -s ./Snakefile --configfile ./config.yaml --cores 1 --unlock --directory /scratch/pawsey0933/cfolland/SVmerge/

#Run snakemake
snakemake --directory /scratch/pawsey0933/cfolland/SVmerge -s ./Snakefile --configfile ./config.yaml --profile svmerge
