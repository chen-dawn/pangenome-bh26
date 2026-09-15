#!/bin/bash
# vg giraffe test of the whole-MHC Minigraph-Cactus graph with one 1000G sample's extracted MHC reads
set -uo pipefail
C=$HOME/hla/cactus/cactus-bin-v3.3.0; export PATH=$C/bin:$HOME/hla/mm/envs/typing/bin:$PATH
O=$HOME/hla/mhcgraph/out; W=$HOME/hla/mhcgraph/giraffe_test; mkdir -p $W && cd $W
vg version | head -1
echo "haplotype samples in GBZ: $(vg paths -M -x $O/MHC.gbz | awk -F'\t' 'NR>1 && $2=="HAPLOTYPE"{print $3}' | sort -u | wc -l)"
echo "haplotype paths: $(vg paths -M -x $O/MHC.gbz | awk -F'\t' 'NR>1 && $2=="HAPLOTYPE"' | wc -l)  reference paths: $(vg paths -M -x $O/MHC.gbz | awk -F'\t' 'NR>1 && $2=="REFERENCE"{print $1}' | tr '\n' ' ')"
R=$HOME/hla/ref1kg/local/GRCh38_full_analysis_set_plus_decoy_hla.fa
S=HG00096; CR=/home/asianhla/data/upload/1000G_MHC/cram/$S.mhc.cram
samtools collate -@${SLURM_CPUS_PER_TASK:-8} -u -O --reference $R $CR tmp | samtools fastq -@${SLURM_CPUS_PER_TASK:-8} -F 0x900 --reference $R -1 r1.fq.gz -2 r2.fq.gz -0 /dev/null -s /dev/null -n - 2>/dev/null
/usr/bin/time -f "giraffe: %e s, %M KB" vg giraffe -t ${SLURM_CPUS_PER_TASK:-8} -Z $O/MHC.gbz -d $O/MHC.dist -m $O/MHC.shortread.withzip.min -z $O/MHC.shortread.zipcodes \
  -f r1.fq.gz -f r2.fq.gz -o gaf > $S.gaf 2> giraffe.err; tail -2 giraffe.err
n=$(($(zcat r1.fq.gz | wc -l)/4*2)); a=$(awk '$6!="*"' $S.gaf | cut -f1 | sort -u | wc -l); q=$(awk '$6!="*" && $12>=20' $S.gaf | cut -f1 | sort -u | wc -l)
echo "reads $n  aligned $a  MAPQ>=20 $q"
echo GIRAFFEDONE
