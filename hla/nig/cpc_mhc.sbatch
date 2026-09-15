#!/bin/bash
#SBATCH -p epyc -A general_analysis
#SBATCH -c 4 --mem=48G -t 24:00:00 -J cpcmhc -o cpcmhc.%j.log
# CPC Phase 1 (Chinese Pangenome Consortium, Minigraph-Cactus graph on CHM13v2; per-sample FASTAs on the cluster are chr16
# only, and CPC.Phase1.CHM13v2-full.odgi carries only the reference paths): cut the MHC from the haplotype walks (W lines)
# of CPC.Phase1.CHM13v2-full.gfa.
#  1. odgi extract of CHM13v2.chr6:28-34 Mb from the odgi graph (node ids = GFA segment ids) -> reference MHC node set
#  2. GFA pass 1: for every CPC sample walk (HIFI*D, RY*) in the chr6 block, keep the steps from the first to the last
#     reference MHC node (+-200 steps); GFA pass 2: sequences of those segments
#  3. write <sample>.hap<N>.all.fa (records <contig>_<walkstart>_<piece>), keep pieces with an alignment of >= 50 kb matching
#     bases (MAPQ >= 20) to the GRCh38 MHC -> fasta/<sample>.hap<N>.mhcctg.fa; extract_mhc.py trims to the MHC later.
# HG00438/HG00621/HG00673 (HPRC, already included) and GRCh38 are skipped.
set -euo pipefail
export PATH=$HOME/hla/mm/envs/hla/bin:$PATH
W=$HOME/hla/cpc; cd $W; mkdir -p fasta tmp
GFA=/home/asianhla/data/CPC/CPC.Phase1.CHM13v2-full.gfa
[ -s mhc_ref.fa ] || samtools faidx $HOME/hla/kg/chr6.fa chr6:28510120-33480577 > mhc_ref.fa
if [ ! -s tmp/ref_mhc_nodes.txt ]; then
  [ -s tmp/mhc.og ] || odgi extract -i /home/asianhla/data/CPC/CPC.Phase1.CHM13v2-full.odgi -r CHM13v2.chr6:28000000-34000000 -c 0 -E -t 4 -P -o tmp/mhc.og
  odgi view -i tmp/mhc.og -g | awk -F'\t' '$1=="P"{n=split($3,a,","); for(i=1;i<=n;i++){sub(/[+-]$/,"",a[i]); print a[i]}}' > tmp/ref_mhc_nodes.txt
fi
echo "reference MHC nodes: $(wc -l < tmp/ref_mhc_nodes.txt)"
python3 - <<'PY'
import re, sys, array, time
GFA = "/home/asianhla/data/CPC/CPC.Phase1.CHM13v2-full.gfa"
keep = re.compile(r"^(HIFI\d+D|RY\d+)$")
ref = set(int(x) for x in open("tmp/ref_mhc_nodes.txt"))
step = re.compile(r"([<>])(\d+)")
PAD = 200
t0 = time.time(); pieces = []; need = set(); block = None
with open(GFA) as f:
    for line in f:
        c = line[0]
        if c == "P":
            block = line.split("\t", 2)[1]
            if block == "CHM13v2.chr7": break
            continue
        if c != "W" or block != "CHM13v2.chr6": continue
        t = line.split("\t", 6)
        if not keep.match(t[1]): continue
        ids = array.array("q", ((int(i) if o == ">" else -int(i)) for o, i in step.findall(t[6])))
        hit = [k for k, v in enumerate(ids) if abs(v) in ref]
        if not hit: continue
        # split at large gaps between reference hits (a walk may pass the MHC once; guard against stray repeats)
        runs, s0, prev = [], hit[0], hit[0]
        for k in hit[1:]:
            if k - prev > 50000: runs.append((s0, prev)); s0 = k
            prev = k
        runs.append((s0, prev))
        for j, (a, b) in enumerate(runs):
            if b - a < 100: continue
            sub = ids[max(0, a - PAD): b + PAD + 1]
            pieces.append((t[1], t[2], t[3], t[4], j, sub)); need.update(abs(v) for v in sub)
print(f"pass 1: {len(pieces)} walk pieces, {len(need)} segments, {time.time()-t0:.0f}s", flush=True)
seq = {}
with open(GFA) as f:
    for line in f:
        if line[0] != "S": continue
        i = line.find("\t", 2); sid = int(line[2:i])
        if sid in need:
            seq[sid] = line[i + 1:].split("\t", 1)[0].rstrip("\n")
            if len(seq) == len(need): break
print(f"pass 2: {len(seq)} sequences, {time.time()-t0:.0f}s", flush=True)
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
out = {}
for s, h, ctg, st, j, sub in pieces:
    fn = f"tmp/{s}.hap{h}.all.fa"
    if fn not in out: out[fn] = open(fn, "w")
    body = "".join(seq[v] if v > 0 else seq[-v].translate(comp)[::-1] for v in sub)
    name = re.sub(r"[^A-Za-z0-9._-]", "_", f"{ctg}_{st}_{j}")
    out[fn].write(f">{name}\n")
    for k in range(0, len(body), 80): out[fn].write(body[k:k + 80] + "\n")
for h in out.values(): h.close()
print(f"{len(out)} CPC haplotype files, {time.time()-t0:.0f}s")
PY
for t in tmp/*.all.fa; do
  [ -e "$t" ] || break
  b=$(basename $t .all.fa); out=fasta/$b.mhcctg.fa
  minimap2 -t 4 -x asm20 --secondary=no mhc_ref.fa $t > tmp/$b.paf 2> tmp/$b.mm2.err
  awk '$10>=50000 && $12>=20{print $1}' tmp/$b.paf | sort -u > tmp/$b.keep
  if [ ! -s tmp/$b.keep ]; then echo "no MHC segment for $b"; continue; fi
  samtools faidx $t; samtools faidx -r tmp/$b.keep $t > $out; samtools faidx $out
  echo "$b $(wc -l < $out.fai) segments $(awk '{s+=$2}END{print s}' $out.fai) bp"
done
echo CPCDONE
