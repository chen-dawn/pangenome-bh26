#!/bin/bash
# Relaunch the workflow from precomputed MHC fastas (skips extraction).
# usage: bash nig/stage2.sh <dir with *.mhc.fa/*.mhc.tsv> [<more dirs>...]
set -euo pipefail
H=$HOME/hla
export PATH=$H/mm/envs/hla/bin:$PATH
MHC=$H/mhc_all; mkdir -p "$MHC"
for d in "$@"; do cp -n "$d"/*.mhc.fa "$d"/*.mhc.tsv "$MHC"/; done
echo "$(ls "$MHC" | grep -c 'mhc.fa$') MHC fastas in $MHC"
cd $H/cwl
python3 scripts/make_inputs.py \
  --dir APR=/home/asianhla/data/upload/APR/assemblies \
  --dir HPRC_r2=/home/asianhla/data/HPRC_r2/fasta \
  --dir JaSaPaGe=/home/asianhla/data/JaSaPaGe/assembly_clean/fasta \
  --reference /home/asianhla/data/HPRC_r2/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz \
  --immuannot-dir $H/Immuannot --immuannot-ref $H/Data-2024Feb02 \
  --mhc-dir "$MHC" --out nig/inputs-nig-stage2.yml 2>&1 | grep -v skipped
cd $H && rm -rf results-stage2 jobstore-results-stage2 && mkdir -p results-stage2
sbatch --cpus-per-task=28 --mem=200G --time=1-00:00:00 --job-name=hla_stage2 \
  --output=results-stage2/slurm-%j.out --error=results-stage2/slurm-%j.err \
  cwl/nig/run_nig.sbatch $H/cwl/nig/inputs-nig-stage2.yml $H/results-stage2
