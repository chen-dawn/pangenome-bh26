#!/usr/bin/env python3
"""Tile assembled MHC haplotype sequences into overlapping error-free pseudo-HiFi reads for FuFiHLA
(which expects diploid long reads). Both haplotypes of an individual go into one file.
usage: tile_reads.py --len 15000 --step 1500 out.fa.gz hap1.mhc.fa [hap2.mhc.fa]"""
import argparse, gzip
ap = argparse.ArgumentParser(); ap.add_argument("--len", type=int, default=15000); ap.add_argument("--step", type=int, default=1500)
ap.add_argument("out"); ap.add_argument("fastas", nargs="+"); a = ap.parse_args()
def records(p):
    name, seq = None, []
    for line in (gzip.open(p, "rt") if p.endswith(".gz") else open(p)):
        if line.startswith(">"):
            if name: yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else: seq.append(line.strip())
    if name: yield name, "".join(seq)
n = 0
with gzip.open(a.out, "wt", compresslevel=1) as out:
    for p in a.fastas:
        for name, s in records(p):
            L = len(s)
            starts = list(range(0, max(1, L - a.len) + 1, a.step))
            if starts[-1] + a.len < L: starts.append(max(0, L - a.len))
            for st in starts:
                out.write(f">r{n}_{name.replace('#','_')}_{st}\n{s[st:st + a.len]}\n"); n += 1
print(n, "pseudo-reads")
