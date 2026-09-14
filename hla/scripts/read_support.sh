#!/bin/bash
# Read support on the assemblies: 1000G high-coverage Illumina reads from the MHC
# region (fetched from the public CRAMs into ../<sample>.mhc.bam) realigned to the
# same individual's own assembled MHC haplotype(s). Produces per-site allele
# counts (sites with an alternative base at >=15%), and an asm20 alignment of
# haplotype 1 against haplotype 2 for each individual.
export PATH=$HOME/hla/mm/envs/hla/bin:$PATH OPENBLAS_NUM_THREADS=1
cd ~/hla/kg/rs
until grep -q ALLDONE ../fetch.log; do sleep 20; done

cat > sites.py <<"PY"
import sys, re
for line in sys.stdin:
    f = line.split("\t")
    if len(f) < 5:
        continue
    d = int(f[3])
    if d < 8:
        continue
    b = re.sub(r"\^.", "", f[4]).replace("$", "")
    # drop indel bases that follow +N/-N
    out, i = [], 0
    while i < len(b):
        c = b[i]
        if c in "+-":
            m = re.match(r"[+-](\d+)", b[i:])
            n = int(m.group(1)); i += len(m.group(0)) + n
            continue
        out.append(c); i += 1
    b = "".join(out)
    ref = sum(1 for c in b if c in ".,")
    alt = {}
    for c in b.upper():
        if c in "ACGT":
            alt[c] = alt.get(c, 0) + 1
    if not alt:
        continue
    a, n = max(alt.items(), key=lambda x: x[1])
    tot = ref + sum(alt.values())
    frac = n / tot if tot else 0
    if frac >= 0.15:
        print(f"{f[0]}\t{f[1]}\t{f[2]}\t{d}\t{ref}\t{a}\t{n}\t{frac:.3f}")
PY

while read s r labs; do
  samtools fastq -N ../$r.mhc.bam 2>/dev/null | gzip -1 > $s.fq.gz
  for lab in $labs; do
    ref=~/hla/mhc_all/${lab}.mhc.fa
    minimap2 -t 8 -ax sr --secondary=no $ref $s.fq.gz 2>/dev/null | samtools sort -@4 -o $s.$lab.bam - && samtools index $s.$lab.bam
    samtools mpileup -q 20 -Q 20 -d 5000 -f $ref $s.$lab.bam 2>/dev/null | python3 sites.py > $s.$lab.sites.tsv
    samtools depth -a $s.$lab.bam | awk '{b=int($2/50000); s[b]+=$3; n[b]++} END{for(k in s) print k*50000, s[k]/n[k]}' | sort -n > $s.$lab.depth50k.tsv
  done
done <<"TAB"
NA18976 NA18976 NA18976_1
NA19909 NA19909 NA19909_1
HG00544 HG00544 HG00544_1
NA18940 NA18940 NA18940_1
NA18952 NA18952 NA18952_1 NA18952_2
NA18952j NA18952 NA18952.JaSaPaGe_1 NA18952.JaSaPaGe_2
HG02717 HG02717 HG02717_1
NA20346 NA20346 NA20346_2
TAB

for s in NA18976 NA19909 HG00544 NA18940 NA18952 NA18952.JaSaPaGe; do
  minimap2 -t 8 -cx asm20 --cs ~/hla/mhc_all/${s}_1.mhc.fa ~/hla/mhc_all/${s}_2.mhc.fa 2>/dev/null > $s.h1h2.paf
done
echo RSDONE
