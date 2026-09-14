#!/usr/bin/env python3
"""Aggregate Immuannot GTFs across haplotypes.

Outputs
  hla_calls.tsv                 one row per annotated gene copy
  hla_calls_matrix.tsv          haplotype x gene matrix of consensus calls
  gene_copy_number.tsv          haplotype x gene copy counts
  mhc_extraction_summary.tsv    concatenated per-haplotype extraction stats
  queries/<GENE>.fa             gene +/- flank from a reference haplotype, query for pgr-query
  genes.txt                     list of genes with a query
  plots/*.png                   allele-frequency bars per gene and cohort, copy-number heatmap,
                                extraction coverage
"""
import argparse
import gzip
import os
import re
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
ATTR = re.compile(r'(\S+) "?([^";]*)"?;')


def attrs(s):
    return dict(ATTR.findall(s))


def read_gtf(path):
    genes, transcripts = {}, {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            key = a.get("gene_id")
            if f[2] == "gene":
                genes[key] = dict(contig=f[0], source=f[1], start=int(f[3]), end=int(f[4]), strand=f[6],
                                  gene=a.get("gene_name"), template_allele=a.get("template_allele", "NA"),
                                  template_distance=a.get("template_distance", "NA"))
            elif f[2] == "transcript":
                transcripts.setdefault(key, dict(consensus=a.get("consensus", a.get("allele", "NA")),
                                                 alleles=a.get("alleles", "NA"),
                                                 template_warning=a.get("template_warning", "NA"),
                                                 cds_distance=a.get("cds_distance", "NA")))
    rows = []
    for k, g in genes.items():
        t = transcripts.get(k, {})
        g.update(consensus=t.get("consensus", "NA"), alleles=t.get("alleles", "NA"),
                 template_warning=t.get("template_warning", "NA"), cds_distance=t.get("cds_distance", "NA"))
        rows.append(g)
    return rows


def fields(allele, n):
    """Trim an allele name to n colon-separated fields; strip Immuannot suffixes."""
    if not isinstance(allele, str) or "*" not in allele:
        return allele
    gene, rest = allele.split("*", 1)
    parts = rest.split(":")
    return gene + "*" + ":".join(parts[:n])


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
    ap.add_argument("--gtf", nargs="+", required=True)
    ap.add_argument("--mhc-fasta", nargs="+", required=True)
    ap.add_argument("--mhc-tsv", nargs="+", required=True)
    ap.add_argument("--sample", nargs="+", required=True)
    ap.add_argument("--haplotype", nargs="+", required=True)
    ap.add_argument("--cohort", nargs="+", required=True)
    ap.add_argument("--genes", default="HLA-A,HLA-B,HLA-C,HLA-E,HLA-F,HLA-G,HLA-DRA,HLA-DRB1,HLA-DRB3,HLA-DRB4,HLA-DRB5,HLA-DQA1,HLA-DQB1,HLA-DPA1,HLA-DPB1,MICA,MICB,TAP1,TAP2,C4A,C4B")
    ap.add_argument("--gene-flank", type=int, default=2000)
    ap.add_argument("--min-graph-seqs", type=int, default=4)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    n = len(a.gtf)
    assert n == len(a.mhc_fasta) == len(a.mhc_tsv) == len(a.sample) == len(a.haplotype) == len(a.cohort)
    os.makedirs(os.path.join(a.outdir, "queries"), exist_ok=True)
    os.makedirs(os.path.join(a.outdir, "plots"), exist_ok=True)

    rows = []
    for i in range(n):
        for g in read_gtf(a.gtf[i]):
            g.update(cohort=a.cohort[i], sample=a.sample[i], haplotype=a.haplotype[i])
            rows.append(g)
    calls = pd.DataFrame(rows)
    if calls.empty:
        raise SystemExit("no gene annotations found")
    calls["novel"] = calls["consensus"].astype(str).str.endswith(":new") | calls["consensus"].astype(str).str.contains("new")
    calls["allele_2field"] = calls["consensus"].map(lambda x: fields(x, 2))
    calls["allele_3field"] = calls["consensus"].map(lambda x: fields(x, 3))
    calls["allele_4field"] = calls["consensus"].map(lambda x: fields(x, 4))

    def novel_class(r):
        if not r["novel"]:
            return "known"
        try:
            return "novel_noncoding" if float(r["cds_distance"]) == 0 else "novel_coding"
        except (TypeError, ValueError):
            return "novel_unknown"
    calls["novel_class"] = calls.apply(novel_class, axis=1)
    calls["hap_id"] = calls["sample"] + "#" + calls["haplotype"]
    cols = ["cohort", "sample", "haplotype", "hap_id", "gene", "contig", "start", "end", "strand", "consensus",
            "allele_2field", "allele_3field", "allele_4field", "novel", "novel_class", "template_allele", "template_distance", "cds_distance",
            "template_warning", "alleles", "source"]
    calls = calls[cols].sort_values(["cohort", "sample", "haplotype", "gene", "start"])
    calls.to_csv(os.path.join(a.outdir, "hla_calls.tsv"), sep="\t", index=False)

    matrix = calls.groupby(["cohort", "hap_id", "gene"])["consensus"].agg(";".join).unstack("gene").fillna("")
    matrix.to_csv(os.path.join(a.outdir, "hla_calls_matrix.tsv"), sep="\t")
    cn = calls.groupby(["cohort", "hap_id", "gene"]).size().unstack("gene").fillna(0).astype(int)
    cn.to_csv(os.path.join(a.outdir, "gene_copy_number.tsv"), sep="\t")

    ext = pd.concat([pd.read_csv(p, sep="\t").assign(cohort=c) for p, c in zip(a.mhc_tsv, a.cohort)], ignore_index=True)
    ext.to_csv(os.path.join(a.outdir, "mhc_extraction_summary.tsv"), sep="\t", index=False)

    # per-gene query sequences for pgr-query: gene +/- flank taken from a reference
    # haplotype (GRCh38 preferred, then CHM13, then the first haplotype carrying the gene)
    wanted = [g for g in a.genes.split(",") if g]
    fastas = {}
    idx = {(a.sample[i], a.haplotype[i]): i for i in range(n)}
    pref = {"GRCh38": 0, "grch38": 0, "chm13": 1, "CHM13": 1}
    written = []
    for gene in wanted:
        sub = calls[calls["gene"] == gene]
        if len(sub) < a.min_graph_seqs:
            continue
        sub = sub.assign(_pref=sub["sample"].map(lambda x: pref.get(x, 2))).sort_values(["_pref", "template_distance"])
        r = next(iter(sub.itertuples()))
        i = idx[(r.sample, r.haplotype)]
        if i not in fastas:
            fastas[i] = read_fasta(a.mhc_fasta[i])
        s_ = fastas[i].get(r.contig)
        if s_ is None:
            continue
        st, en = max(0, r.start - 1 - a.gene_flank), min(len(s_), r.end + a.gene_flank)
        seq = s_[st:en]
        if r.strand == "-":
            seq = seq.translate(COMP)[::-1]
        with open(os.path.join(a.outdir, "queries", f"{gene}.query.fa"), "w") as out:
            out.write(f">{gene} source={r.sample}#{r.haplotype}#{r.contig}:{st + 1}-{en}({r.strand}) template={r.template_allele}\n")
            for j in range(0, len(seq), 80):
                out.write(seq[j:j + 80] + "\n")
        written.append(gene)
    with open(os.path.join(a.outdir, "genes.txt"), "w") as fh:
        fh.write("\n".join(written) + ("\n" if written else ""))

    # plots: allele frequencies per gene by cohort (2-field)
    cohorts = sorted(calls["cohort"].unique())
    for gene in wanted:
        sub = calls[(calls["gene"] == gene) & (calls["cohort"] != "REF")]
        if sub.empty:
            continue
        freq = sub.groupby(["cohort", "allele_2field"]).size().unstack("cohort").fillna(0)
        freq = freq / freq.sum()
        top = freq.max(axis=1).sort_values(ascending=False).head(25).index
        freq = freq.loc[top]
        fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(top) + 2), 4.5))
        x = np.arange(len(top))
        w = 0.8 / max(1, len(freq.columns))
        for j, c in enumerate(freq.columns):
            ax.bar(x + j * w, freq[c].values, w, label=f"{c} (n={int((sub['cohort'] == c).sum())})")
        ax.set_xticks(x + w * (len(freq.columns) - 1) / 2)
        ax.set_xticklabels([t.split("*", 1)[1] if "*" in t else t for t in top], rotation=90, fontsize=8)
        ax.set_ylabel("haplotype frequency")
        ax.set_title(f"{gene} two-field alleles (Immuannot consensus)")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "plots", f"allele_freq_{gene}.png"), dpi=130)
        plt.close(fig)

    # copy-number heatmap (mean copies per haplotype by cohort)
    cn2 = calls[calls["cohort"] != "REF"].groupby(["cohort", "hap_id", "gene"]).size().unstack("gene").fillna(0)
    haps_per_cohort = calls[calls["cohort"] != "REF"].groupby("cohort")["hap_id"].nunique()
    mean_cn = cn2.groupby(level="cohort").sum().div(haps_per_cohort, axis=0)
    fig, ax = plt.subplots(figsize=(max(8, 0.3 * mean_cn.shape[1]), 1 + 0.5 * mean_cn.shape[0] + 2))
    im = ax.imshow(mean_cn.values, aspect="auto", cmap="viridis", vmin=0, vmax=max(1.0, float(mean_cn.values.max())))
    ax.set_yticks(range(mean_cn.shape[0]))
    ax.set_yticklabels([f"{c} (n={haps_per_cohort[c]})" for c in mean_cn.index])
    ax.set_xticks(range(mean_cn.shape[1]))
    ax.set_xticklabels(mean_cn.columns, rotation=90, fontsize=7)
    fig.colorbar(im, ax=ax, label="mean copies per haplotype")
    ax.set_title("HLA/KIR/C4 gene presence per haplotype (Immuannot)")
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "plots", "gene_copy_number.png"), dpi=130)
    plt.close(fig)

    # novel-allele rate per gene
    nc = calls[calls["cohort"] != "REF"]
    nov = nc.groupby(["gene", "cohort"])["novel"].mean().unstack("cohort").fillna(0)
    nov = nov.loc[[g for g in wanted if g in nov.index]]
    if not nov.empty:
        fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(nov) + 2), 4))
        nov.plot.bar(ax=ax, width=0.8)
        ax.set_ylabel("fraction of gene copies with novel allele")
        ax.set_title("Alleles absent from IPD-IMGT/HLA 3.55 (Immuannot 'new')")
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "plots", "novel_allele_rate.png"), dpi=130)
        plt.close(fig)
        cls = nc[nc["novel"]].groupby(["gene", "novel_class"]).size().unstack("novel_class").fillna(0)
        cls = cls.loc[[g for g in wanted if g in cls.index]]
        fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(cls) + 2), 4))
        cls.plot.bar(ax=ax, stacked=True, width=0.8)
        ax.set_ylabel("novel gene copies")
        ax.set_title("Novel alleles: coding (CDS differs) vs noncoding-only (cf. Ito-Naito et al. 2026)")
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "plots", "novel_allele_class.png"), dpi=130)
        plt.close(fig)

    # extraction coverage
    cov = ext.groupby(["cohort", "sample", "haplotype"])["mhc_covered_frac"].sum().reset_index()
    fig, ax = plt.subplots(figsize=(7, 4))
    data = [cov.loc[cov["cohort"] == c, "mhc_covered_frac"].values for c in cohorts]
    ax.boxplot(data, tick_labels=cohorts)
    ax.set_ylabel("fraction of GRCh38 MHC covered by alignments")
    ax.set_title("MHC region extraction per haplotype")
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "plots", "mhc_extraction_coverage.png"), dpi=130)
    plt.close(fig)

    # diversity: distinct 2-field alleles per gene per cohort
    div = calls[calls["cohort"] != "REF"].groupby(["gene", "cohort"])["allele_2field"].nunique().unstack("cohort").fillna(0)
    div = div.loc[[g for g in wanted if g in div.index]]
    if not div.empty:
        fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(div) + 2), 4))
        div.plot.bar(ax=ax, width=0.8)
        ax.set_ylabel("distinct two-field alleles")
        ax.set_title("Allelic diversity per gene and cohort")
        fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "plots", "allele_diversity.png"), dpi=130)
        plt.close(fig)


if __name__ == "__main__":
    main()
