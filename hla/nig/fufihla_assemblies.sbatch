#!/bin/bash
#SBATCH -p epyc -A general_analysis
#SBATCH -c 4 --mem=16G -t 2:00:00 -J fufihla
#SBATCH -o /home/asianhla/data/upload/HLA/fufihla/logs/%A_%a.log
# FuFiHLA (Hu et al., Bioinformatics 2026, btag231; full-field typing of HLA-A/B/C/DRB1/DQA1/DQB1 from long reads) on the
# assembled MHC haplotypes: both haplotypes of an individual are tiled into error-free 15 kb pseudo-HiFi reads every
# 1.5 kb (10x each, scripts/tile_reads.py) and typed together as one diploid sample. Checked on HG002: 12/12 calls equal Immuannot.
# individuals.tsv: sample <tab> hap1 MHC fasta <tab> hap2 MHC fasta
# usage: sbatch --array=1-N%40 fufihla_assemblies.sbatch
set -euo pipefail
O=/home/asianhla/data/upload/HLA/fufihla
export PATH=$HOME/hla/mm/envs/typing/bin:$PATH
read -r s f1 f2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $O/individuals.tsv)
[ -s $O/calls/$s.out ] && { echo "$s already typed"; exit 0; }
T=${SLURM_TMPDIR:-/tmp}/fufi_$s; rm -rf $T; mkdir -p $T $O/calls $O/consensus
python3 $HOME/hla/cwl/scripts/tile_reads.py $T/reads.fa.gz $f1 $f2
cd $T && fufihla --fa reads.fa.gz --out $s --refdir $HOME/hla/typing/ref_data > $s.log 2>&1
cp $s.out $O/calls/$s.out; cp $s.log $O/calls/$s.log; cp -r $s/consensus $O/consensus/$s
cd / && rm -rf $T
echo "$s done $(date)"
