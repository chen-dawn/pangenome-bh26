#!/usr/bin/env python3
"""Post-process pgr-query output: filter hits, rename records to PanSN names.

pgr-query writes <prefix>.000.hit (one row per target hit) and <prefix>.000.fa
with records named <db_stem>::<contig>_<start>_<end>_<orientation>. Sequences
with orientation 1 are already reverse-complemented by pgr-query, so every
output record is in query orientation.

Two naming modes:
  --sample S --haplotype H   whole-assembly extraction: names S#H#contig:start-end[_rc]
  --gene G                   fetch from PanSN-named database: names sample#hap#G[_n]
"""
import argparse
import sys
from collections import Counter


def read_hits(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            # observed column order (pgr-tk 0.5.1): idx, q_name, q_bgn, q_end, q_len,
            # anchors, source, ctg, ctg_bgn, ctg_end, orientation, out_seq_name
            rows.append(dict(q_name=f[1], q_bgn=int(f[2]), q_end=int(f[3]), q_len=int(f[4]),
                             anchors=int(f[5]), source=f[6], ctg=f[7], ctg_bgn=int(f[8]),
                             ctg_end=int(f[9]), orientation=int(f[10]), name=f[11]))
    return rows


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hit", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--sample")
    ap.add_argument("--haplotype")
    ap.add_argument("--gene")
    ap.add_argument("--min-piece", type=int, default=20000, help="minimum target span to keep")
    ap.add_argument("--min-anchors", type=int, default=20)
    ap.add_argument("--min-query-frac", type=float, default=0.0, help="minimum fraction of the query covered")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-tsv", required=True)
    a = ap.parse_args()
    if not a.gene and not (a.sample and a.haplotype):
        sys.exit("need --gene or --sample/--haplotype")

    hits = read_hits(a.hit)
    seqs = read_fasta(a.fasta)
    keep = []
    for h in hits:
        span = h["ctg_end"] - h["ctg_bgn"]
        qfrac = (h["q_end"] - h["q_bgn"]) / h["q_len"] if h["q_len"] else 0
        if span < a.min_piece or h["anchors"] < a.min_anchors or qfrac < a.min_query_frac:
            continue
        keep.append(h)
    keep.sort(key=lambda h: -(h["q_end"] - h["q_bgn"]))

    copies = Counter()
    q_len = hits[0]["q_len"] if hits else 0
    covered = sum(h["q_end"] - h["q_bgn"] for h in keep)
    with open(a.out_fasta, "w") as fa, open(a.out_tsv, "w") as tsv:
        tsv.write("sample\thaplotype\tcontig\tcontig_len\tstart\tend\tstrand\tpiece_len\tmhc_covered_bp\tmhc_covered_frac\tn_pieces\tanchors\tquery\n")
        for h in keep:
            seq = seqs.get(h["name"])
            if seq is None:
                print(f"missing sequence {h['name']}", file=sys.stderr)
                continue
            strand = "-" if h["orientation"] == 1 else "+"
            if a.gene:
                parts = h["ctg"].split("#")
                sample, hap = (parts[0], parts[1]) if len(parts) >= 3 else (h["ctg"], "0")
                copies[(sample, hap)] += 1
                suffix = f"_{copies[(sample, hap)]}" if copies[(sample, hap)] > 1 else ""
                name = f"{sample}#{hap}#{a.gene}{suffix}"
            else:
                sample, hap = a.sample, a.haplotype
                ctg = h["ctg"]
                # assemblies that are already PanSN-named (HPRC r2, references) would
                # otherwise get a second sample#hap# prefix; the sample may carry a
                # ".cohort" suffix added to disambiguate duplicates across cohorts
                for base in (sample, sample.rsplit(".", 1)[0]):
                    pre = f"{base}#{hap}#".lower()
                    if ctg.lower().startswith(pre):
                        ctg = ctg[len(pre):]
                        break
                name = f"{sample}#{hap}#{ctg}:{h['ctg_bgn'] + 1}-{h['ctg_end']}" + ("_rc" if strand == "-" else "")
            fa.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i + 80] + "\n")
            tsv.write(f"{sample}\t{hap}\t{h['ctg']}\tNA\t{h['ctg_bgn'] + 1}\t{h['ctg_end']}\t{strand}\t{h['ctg_end'] - h['ctg_bgn']}\t{h['q_end'] - h['q_bgn']}\t{(h['q_end'] - h['q_bgn']) / q_len if q_len else 0:.4f}\t{len(keep)}\t{h['anchors']}\t{h['q_name']}\n")
        if not keep and not a.gene:
            tsv.write(f"{a.sample}\t{a.haplotype}\tNA\tNA\t0\t0\tNA\t0\t0\t0.0000\t0\t0\tNA\n")
    print(f"kept {len(keep)}/{len(hits)} hits, query covered {covered / q_len if q_len else 0:.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
