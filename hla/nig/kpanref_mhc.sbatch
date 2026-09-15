#!/bin/bash
#SBATCH --partition=asianhla-c32 --account=asianhla-group
#SBATCH -c 24 --mem=150G -t 24:00:00 -J kpmhc -o kpmhc.%j.log
# K-PanRef (Korean pangenome, Minigraph-Cactus GBZ, no public assembly FASTAs): pull every haplotype path of each
# sample from the graph and keep the path segments that align to the GRCh38 MHC (one alignment with >= 50 kb matching bases, MAPQ >= 20, minimap2 asm20).
# Only KOREF1 has chromosome-named contigs in the graph; the others are named Unplaced_h?tg*, so we cannot select chr6 by name.
# Output fasta/<sample>.hap<N>.mhcctg.fa, records <locus>_<offset> (extract_mhc.py adds the sample#hap# prefix); KOREF1-G-TTAGGA -> KOREF1, KOREF2-G -> KOREF2.
set -euo pipefail
export PATH=$HOME/hla/mm/envs/hla/bin:$PATH
W=$HOME/hla/kpanref; cd $W; mkdir -p fasta tmp
G=/home/asianhla/data/upload/KPanRef/KPanRef.gbz
[ -s paths.meta.tsv ] || vg paths -x $G -M > paths.meta.tsv          # path metadata, ~9 GB RAM
[ -s mhc_ref.fa ] || samtools faidx $HOME/hla/kg/chr6.fa chr6:28510120-33480577 > mhc_ref.fa
awk -F'\t' '$2=="HAPLOTYPE"{print $3"\t"$4}' paths.meta.tsv | sort -u > haps.tsv
one() {
  s=$1; h=$2; short=$s; [ "$s" = KOREF1-G-TTAGGA ] && short=KOREF1; [ "$s" = KOREF2-G ] && short=KOREF2
  out=fasta/$short.hap$h.mhcctg.fa; [ -s $out ] && return 0
  t=tmp/$short.$h.fa
  awk -F'\t' -v s=$s -v h=$h '$2=="HAPLOTYPE" && $3==s && $4==h{print $1}' paths.meta.tsv > tmp/$short.$h.paths
  vg paths -x $G -F -p tmp/$short.$h.paths > $t
  minimap2 -t 6 -x asm20 --secondary=no mhc_ref.fa $t > tmp/$short.$h.paf 2> tmp/$short.$h.mm2.err
  # a genuine MHC segment has an alignment with >= 50 kb matching bases (the KOREF1 chr6 segment: one 5.1 Mb block,
  # 5.37 M matches); interspersed repeats give long chained blocks with ~80 matching bases (column 10), so filter on matches
  awk '$10>=50000 && $12>=20{print $1}' tmp/$short.$h.paf | sort -u > tmp/$short.$h.keep
  [ -s tmp/$short.$h.keep ] || { echo "no MHC segment for $short#$h" >&2; return 1; }
  samtools faidx $t
  samtools faidx -r tmp/$short.$h.keep $t | awk '/^>/{split(substr($1,2),a,"#"); print ">"a[3]"_"a[4]; next}{print}' > $out.tmp
  mv $out.tmp $out; samtools faidx $out; rm -f $t $t.fai
  echo "$short#$h $(wc -l < $out.fai) segments $(awk '{s+=$2}END{print s}' $out.fai) bp"
}
export -f one; export G
rm -f fasta/*.chr6.fa fasta/*.chr6.fa.fai
xargs -P 4 -L 1 -r bash -c 'one $0 $1' < haps.tsv
echo MHCDONE
