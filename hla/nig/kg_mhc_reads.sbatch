#!/bin/bash
#SBATCH -c 4 --mem=12G -t 4:00:00 -J kgmhc
#SBATCH -o /home/asianhla/data/upload/1000G_MHC/logs/%A_%a.log
# 1000 Genomes high-coverage (NYGC, 2504 unrelated samples; CRAMs aligned to GRCh38_full_analysis_set_plus_decoy_hla):
# per sample, keep every read aligned to the MHC (GRCh38 chr6:28510120-33480577, the region used for assembly
# extraction), to the chr6 alt haplotypes (incl. the MHC alts GL000250-256) and to the 525 HLA-* allele contigs,
# as CRAM against the same reference; then type HLA from those reads with T1K (IPD-IMGT/HLA index).
# usage: sbatch -p <partition> -A <account> --array=1-2504%20 kg_mhc_reads.sbatch
set -euo pipefail
O=/home/asianhla/data/upload/1000G_MHC
D=/usr/local/shared_data/public-human-genomes/GRCh38/1000Genomes/CRAM
R=$HOME/hla/ref1kg/GRCh38_full_analysis_set_plus_decoy_hla.fa    # symlink + own .fai (shared folder is read-only)
IDX=$HOME/hla/typing/t1k/hlaidx
export PATH=$HOME/hla/mm/envs/typing/bin:$PATH
s=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $O/samples.txt)
mkdir -p $O/cram $O/t1k
if [ ! -s $O/cram/$s.mhc.cram.crai ]; then
  samtools view -@3 -C -T $R -o $O/cram/$s.mhc.cram.tmp $D/$s/$s.cram $(cat $O/regions.txt)
  mv $O/cram/$s.mhc.cram.tmp $O/cram/$s.mhc.cram; samtools index $O/cram/$s.mhc.cram
fi
if [ ! -s $O/t1k/$s/${s}_genotype.tsv ]; then
  T=${SLURM_TMPDIR:-/tmp}/kgmhc_$s; mkdir -p $T $O/t1k/$s
  samtools collate -@3 -u -O --reference $R $O/cram/$s.mhc.cram $T/c | samtools fastq -@3 -F 0x900 --reference $R \
    -1 $T/r1.fq.gz -2 $T/r2.fq.gz -0 /dev/null -s /dev/null -n -
  run-t1k -1 $T/r1.fq.gz -2 $T/r2.fq.gz --preset hla -f $IDX/hlaidx_dna_seq.fa -t 4 -o $s --od $O/t1k/$s
  rm -rf $T
fi
echo "$s done $(date)"
