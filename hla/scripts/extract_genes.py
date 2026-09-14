#!/usr/bin/env python3
"""Cut every annotated copy of each requested gene out of the concatenated MHC
haplotype fasta, using the Immuannot coordinates collected in hla_calls.tsv.

Writes genes/<GENE>.fa with records named sample#haplotype#GENE (a second copy
on the same haplotype gets the suffix _2, ...) in gene orientation, gene +/-
flank clipped to the contig, plus genes/<GENE>.regions.tsv listing the cut
regions. Genes with fewer than --min-seqs annotated copies are skipped.
"""
import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True, help="hla_calls.tsv from aggregate_immuannot.py")
    ap.add_argument("--fasta", required=True, help="concatenated MHC haplotypes (PanSN names)")
    ap.add_argument("--genes", required=True, help="comma-separated gene names")
    ap.add_argument("--flank", type=int, default=2000)
    ap.add_argument("--min-seqs", type=int, default=4)
    ap.add_argument("--outdir", default="genes")
    a = ap.parse_args()

    fai = a.fasta + ".fai"
    if not os.path.exists(fai):
        subprocess.run(["samtools", "faidx", a.fasta], check=True)
    ctg_len = {}
    with open(fai) as fh:
        for line in fh:
            f = line.split("\t")
            ctg_len[f[0]] = int(f[1])

    wanted = [g for g in a.genes.split(",") if g]
    rows = defaultdict(list)
    with open(a.calls) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene"] in wanted and r["contig"] in ctg_len:
                rows[r["gene"]].append(r)

    os.makedirs(a.outdir, exist_ok=True)
    written = []
    for gene in wanted:
        recs = rows.get(gene, [])
        if len(recs) < a.min_seqs:
            print(f"{gene}: {len(recs)} copies, skipped", file=sys.stderr)
            continue
        recs.sort(key=lambda r: (r["cohort"], r["sample"], r["haplotype"], int(r["start"])))
        seen = defaultdict(int)
        plan = []  # (name, region, strand, row)
        for r in recs:
            hap = f"{r['sample']}#{r['haplotype']}"
            seen[hap] += 1
            name = f"{hap}#{gene}" + (f"_{seen[hap]}" if seen[hap] > 1 else "")
            s = max(1, int(r["start"]) - a.flank)
            e = min(ctg_len[r["contig"]], int(r["end"]) + a.flank)
            plan.append((name, f"{{{r['contig']}}}:{s}-{e}", r["strand"], r))  # braces: contig names contain ':'

        seqs = {}
        for strand in ("+", "-"):
            regs = [p[1] for p in plan if p[2] == strand]
            if not regs:
                continue
            rf = os.path.join(a.outdir, f"{gene}.{'plus' if strand == '+' else 'minus'}.regions")
            with open(rf, "w") as fh:
                fh.write("\n".join(regs) + "\n")
            cmd = ["samtools", "faidx", "-r", rf] + (["-i"] if strand == "-" else []) + [a.fasta]
            out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
            name = None
            for line in out.splitlines():
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    if name.endswith("/rc"):
                        name = name[:-3]
                    seqs[name] = []
                else:
                    seqs[name].append(line)
            os.remove(rf)
        with open(os.path.join(a.outdir, f"{gene}.fa"), "w") as fa, \
                open(os.path.join(a.outdir, f"{gene}.regions.tsv"), "w") as tsv:
            tsv.write("name\tcohort\tsample\thaplotype\tregion\tstrand\tconsensus\n")
            n = 0
            for name, reg, strand, r in plan:
                if reg not in seqs:
                    print(f"{gene}: region {reg} missing from samtools output", file=sys.stderr)
                    continue
                fa.write(f">{name} {reg.replace('{', '').replace('}', '')}({strand}) allele={r['consensus']}\n")
                fa.write("\n".join(seqs[reg]) + "\n")
                tsv.write(f"{name}\t{r['cohort']}\t{r['sample']}\t{r['haplotype']}\t{reg.replace('{', '').replace('}', '')}\t{strand}\t{r['consensus']}\n")
                n += 1
        written.append((gene, n))
    for gene, n in written:
        print(f"{gene}: {n} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
