#!/bin/bash
#SBATCH -p medium -A general_analysis
#SBATCH -c 32 --mem=240G -t 3-00:00:00 -J mhc_mc -o mhc_mc.%j.log
# Whole-MHC Minigraph-Cactus pangenome of every haplotype (assemblies: APR, HPRC r2, JaSaPaGe; graph-derived: K-PanRef,
# CPC; references GRCh38 and CHM13), for short-read alignment (vg giraffe) of the 1000 Genomes MHC reads.
# Input: haps.tsv (sample, haplotype, MHC fasta from the workflow's extract_mhc step; GRCh38 = GRCh38#0#chr6:28510120-33480577).
# Contigs are renamed to ctg<N> (PanSN '#' and ':' in the original names clash with Cactus/vg path naming); the mapping is
# in contig_names.tsv. Outputs in out/: MHC.gfa.gz / MHC.full.gbz (unclipped), MHC.gbz + MHC.dist + MHC.min (clipped,
# giraffe indexes; GRCh38 reference), MHC.raw.vcf.gz / MHC.vcf.gz, MHC.og (odgi), stats.
set -euo pipefail
W=$HOME/hla/mhcgraph; C=$HOME/hla/cactus/cactus-bin-v3.3.0
cd $W; mkdir -p fa out
. $C/venv-cactus-v3.3.0/bin/activate
export PATH=$C/bin:$PATH PYTHONPATH=$C/lib:${PYTHONPATH:-} TMPDIR=${SLURM_TMPDIR:-/tmp}
if [ ! -s seqfile.txt ]; then
  : > contig_names.tsv; : > seqfile.txt
  while IFS=$'\t' read -r s h f; do
    name=$s; [ "$h" != 0 ] && name=$s.$h
    awk -v s="$s" -v h="$h" -v m=contig_names.tsv '/^>/{n++; print s"\t"h"\tctg"n"\t"substr($1,2) >> m; print ">ctg"n; next}{print}' "$f" > fa/$name.fa
    echo -e "$name\t$W/fa/$name.fa" >> seqfile.txt
  done < haps.tsv
fi
echo "haplotypes: $(wc -l < seqfile.txt)"
cactus-pangenome $W/js $W/seqfile.txt --outDir $W/out --outName MHC --reference GRCh38 CHM13 \
  --batchSystem single_machine --maxCores 32 --maxMemory 230G --workDir $TMPDIR \
  --mapCores 16 --consCores 32 --indexCores 32 \
  --gfa clip full --gbz clip full --giraffe clip --vcf --odgi full --viz --logFile $W/out/cactus.log
vg stats -lz out/MHC.gbz > out/MHC.stats.txt 2>&1 || true
echo MCDONE $(date)
