#!/bin/bash
# Stage 3 on the NIG supercomputer: rerun aggregation, per-gene extraction,
# bundle decomposition and graphs from precomputed MHC fastas and Immuannot GTFs
# (skips extraction and annotation). Usage:
#   bash cwl/nig/stage3.sh <mhc_dir> <gtf_dir> [outname]
# <mhc_dir> holds <sample>_<hap>.mhc.fa/.mhc.tsv, <gtf_dir> holds <sample>_<hap>.gtf.gz.
set -euo pipefail
H=$HOME/hla
MHC=${1:?mhc dir}; GTF=${2:?gtf dir}; OUT=${3:-results-stage3}
export PATH=$H/mm/envs/hla/bin:$PATH
cd $H/cwl
python3 scripts/make_inputs.py \
  --dir APR=/home/asianhla/data/upload/APR/assemblies \
  --dir HPRC_r2=/home/asianhla/data/HPRC_r2/fasta \
  --dir JaSaPaGe=/home/asianhla/data/JaSaPaGe/assembly_clean/fasta \
  --reference /home/asianhla/data/HPRC_r2/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz \
  --immuannot-dir $H/Immuannot --immuannot-ref $H/Data-2024Feb02 \
  --mhc-dir "$MHC" --gtf-dir "$GTF" --out nig/inputs-nig-$OUT.yml 2>&1 | grep -v skipped
cd $H && rm -rf $OUT jobstore-$OUT && mkdir -p $OUT
sbatch --partition=asianhla-c32 --account=asianhla-group \
  --cpus-per-task=28 --mem=200G --time=1-00:00:00 --job-name=hla_$OUT \
  --output=$OUT/slurm-%j.out --error=$OUT/slurm-%j.err \
  cwl/nig/run_nig.sbatch $H/cwl/nig/inputs-nig-$OUT.yml $H/$OUT
