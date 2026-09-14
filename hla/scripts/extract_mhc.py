#!/usr/bin/env python3
"""Extract the MHC region from a haplotype assembly.

Aligns the assembly (query) to a single reference chromosome (target) with
minimap2 (asm preset), projects the reference MHC interval onto the assembly
contigs through the alignments, and writes the covering contig segments in
reference orientation using PanSN-style names ``sample#hap#contig:start-end``.
"""
import argparse
import gzip
import subprocess
import sys
from collections import defaultdict

COMP = str.maketrans("ACGTacgtNnRYKMSWBDHVrykmswbdhv", "TGCAtgcaNnYRMKSWVHDByrmkswvhdb")


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_region(region):
    chrom, span = region.rsplit(":", 1)
    start, end = span.split("-")
    return chrom, int(start.replace(",", "")), int(end.replace(",", ""))


def run_minimap2(ref, asm, threads, preset, paf):
    cmd = ["minimap2", "-t", str(threads), "-x", preset, "--secondary=no", ref, asm]
    with open(paf, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)


def project(qs, qe, strand, ts, te, cs, ce):
    """Linearly project clipped target interval [cs,ce] onto query coords."""
    tl = te - ts
    f1 = (cs - ts) / tl
    f2 = (ce - ts) / tl
    ql = qe - qs
    if strand == "+":
        return int(qs + f1 * ql), int(qs + f2 * ql)
    return int(qs + (1 - f2) * ql), int(qs + (1 - f1) * ql)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--reference", required=True, help="single-chromosome fasta or .mmi")
    ap.add_argument("--region", required=True, help="chrom:start-end (1-based, inclusive)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--haplotype", required=True)
    ap.add_argument("--flank", type=int, default=100000)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-aln", type=int, default=20000)
    ap.add_argument("--merge-gap", type=int, default=500000)
    ap.add_argument("--min-piece", type=int, default=20000)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--preset", default="asm20")
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-paf", default="alignment.paf")
    a = ap.parse_args()

    chrom, rs, re_ = parse_region(a.region)
    rs0 = rs - 1  # 0-based half open
    ws, we = max(0, rs0 - a.flank), re_ + a.flank

    run_minimap2(a.reference, a.assembly, a.threads, a.preset, a.out_paf)

    per_contig = defaultdict(list)  # contig -> list of (qs, qe, strand, weight, covered_target)
    qlen = {}
    with open(a.out_paf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            qname, ql, qs, qe, strand, tname, tl, ts, te, nmatch, alen, mapq = f[:12]
            if tname != chrom:
                continue
            qs, qe, ts, te, alen, mapq = int(qs), int(qe), int(ts), int(te), int(alen), int(mapq)
            if mapq < a.min_mapq or alen < a.min_aln:
                continue
            cs, ce = max(ts, ws), min(te, we)
            if ce - cs <= 0:
                continue
            pqs, pqe = project(qs, qe, strand, ts, te, cs, ce)
            qlen[qname] = int(ql)
            core = max(0, min(te, re_) - max(ts, rs0))
            per_contig[qname].append((pqs, pqe, strand, ce - cs, core))

    pieces = []
    for contig, ivs in per_contig.items():
        ivs.sort()
        merged = []
        for pqs, pqe, strand, w, core in ivs:
            if merged and pqs - merged[-1][1] <= a.merge_gap:
                m = merged[-1]
                m[1] = max(m[1], pqe)
                m[2][strand] += w
                m[3] += core
            else:
                merged.append([pqs, pqe, {"+": 0, "-": 0}, core])
                merged[-1][2][strand] += w
        for pqs, pqe, sw, core in merged:
            if pqe - pqs < a.min_piece:
                continue
            strand = "+" if sw["+"] >= sw["-"] else "-"
            pieces.append((contig, pqs, pqe, strand, core))

    wanted = {p[0] for p in pieces}
    seqs = {}
    if wanted:
        name, buf = None, []
        with open_any(a.assembly) as fh:
            for line in fh:
                if line.startswith(">"):
                    if name in wanted:
                        seqs[name] = "".join(buf)
                    name = line[1:].split()[0]
                    buf = []
                elif name in wanted:
                    buf.append(line.strip())
            if name in wanted:
                seqs[name] = "".join(buf)

    mhc_len = re_ - rs0
    total_core = sum(p[4] for p in pieces)
    with open(a.out_fasta, "w") as fa, open(a.out_tsv, "w") as tsv:
        tsv.write("sample\thaplotype\tcontig\tcontig_len\tstart\tend\tstrand\tpiece_len\tmhc_covered_bp\tmhc_covered_frac\tn_pieces\n")
        for contig, pqs, pqe, strand, core in sorted(pieces, key=lambda p: -p[4]):
            seq = seqs[contig][pqs:pqe]
            if strand == "-":
                seq = revcomp(seq)
            hdr = f"{a.sample}#{a.haplotype}#{contig}:{pqs + 1}-{pqe}"
            if strand == "-":
                hdr += "_rc"
            fa.write(f">{hdr}\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i + 80] + "\n")
            tsv.write(f"{a.sample}\t{a.haplotype}\t{contig}\t{qlen[contig]}\t{pqs + 1}\t{pqe}\t{strand}\t{pqe - pqs}\t{core}\t{core / mhc_len:.4f}\t{len(pieces)}\n")
        if not pieces:
            tsv.write(f"{a.sample}\t{a.haplotype}\tNA\t0\t0\t0\tNA\t0\t0\t0.0000\t0\n")
    print(f"{a.sample}#{a.haplotype}: {len(pieces)} piece(s), MHC covered {total_core / mhc_len:.3f}", file=sys.stderr)


if __name__ == "__main__":
    main()
