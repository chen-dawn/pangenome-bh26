#!/usr/bin/env python3
"""Concordance of HLA typing methods, per individual and gene (HLA-A/B/C/DRB1/DQA1/DQB1):
 - Immuannot on the assembled MHC haplotypes (hla_calls.tsv, column consensus; ':new' suffix dropped)
 - FuFiHLA on tiled pseudo-reads of both assembled haplotypes (fufihla/calls/<sample>.out)
 - T1K on 1000G high-coverage Illumina reads of the MHC (1000G_MHC/t1k/<sample>/<sample>_genotype.tsv), for 1000G
   individuals that were also assembled
Genotypes are compared as unordered allele pairs truncated to 1..4 fields (a pair matches if either pairing matches).
usage: typing_concordance.py --calls hla_calls.tsv --fufihla DIR --t1k DIR --out-prefix data/typing_concordance"""
import argparse, glob, os, re, subprocess, tempfile
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from cohorts import SHORT, present

GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1"]
ap = argparse.ArgumentParser()
ap.add_argument("--calls", default="hla_calls.tsv"); ap.add_argument("--fufihla", required=True)
ap.add_argument("--gene-fasta-dir", help="directory with <GENE>.fa of all annotated assembled gene copies (records sample#hap#GENE)")
ap.add_argument("--t1k"); ap.add_argument("--published", help="1000G 20140702_hla_diversity.txt (Gourraud et al. 2014)")
ap.add_argument("--hprc-truth", help="Lai et al. 2024 HPRC 4-field labels (ground_truth/HPRC_2024/*.xlsx)"); ap.add_argument("--out-prefix", default="data/typing_concordance"); ap.add_argument("--fig", default="figures/fig12_typing_concordance.png")
a = ap.parse_args()

def norm(x):
    if not isinstance(x, str) or "*" not in x: return None
    x = x.strip().split()[0].replace(":new", "").rstrip("#*")      # footnote markers in published tables
    x = re.sub(r"(:00)+$", "", x)                                   # Immuannot pads novel calls, e.g. DRB1*04:92:00
    x = re.sub(r"[A-Z]$", "", x) if re.search(r":\d+[A-Z]$", x) else x   # expression suffix N/L/Q
    return x if x.startswith("HLA-") else "HLA-" + x
def trunc(x, k): return None if x is None else x.split("*")[0] + "*" + ":".join(x.split("*")[1].split(":")[:k])
def pair_match(p, q, k):
    if len(p) != 2 or len(q) != 2 or None in p or None in q: return np.nan
    p = [trunc(x, k) for x in p]; q = [trunc(x, k) for x in q]
    return float((p[0] == q[0] and p[1] == q[1]) or (p[0] == q[1] and p[1] == q[0]))
def depth(p):
    return min(len(x.split("*")[1].split(":")) for x in p) if len(p) == 2 and None not in p else 0

d = pd.read_csv(a.calls, sep="\t", low_memory=False)
d = d[d.gene.isin(GENES) & (d.cohort != "REF")]
imm = {}; imm_novel = {}
for (s, g), x in d.groupby(["sample", "gene"]):
    by_hap = x.sort_values("start").groupby("haplotype").consensus.first()
    imm[(s, g)] = [norm(by_hap.get(h)) for h in sorted(by_hap.index)]
    # Immuannot marks sequences absent from its IPD-IMGT/HLA release with ':new' or reports only the allele group
    imm_novel[(s, g)] = any(isinstance(v, str) and (v.endswith(":new") or ":" not in v) for v in by_hap.values)
cohort = d.groupby("sample").cohort.first().to_dict()

fufi = {}; fufi_inexact = {}
for f in glob.glob(os.path.join(a.fufihla, "*.out")):
    s = os.path.basename(f)[:-4]; al = {}
    for line in open(f):
        t = line.rstrip("\n").split("\t")
        if not t[0].startswith("HLA-") or "*" not in t[0]: continue
        g = t[0].split("*")[0]; al.setdefault(g, []).append(norm(t[0]))
        # exact reconstruction: cs:Z::<length>; anything else = consensus differs from the named allele (novel)
        if not re.fullmatch(r"cs:Z::\d+", t[-1]): fufi_inexact[(s, g)] = True
    for g, v in al.items(): fufi[(s, g)] = v if len(v) == 2 else (v * 2 if len(v) == 1 else v[:2])

t1k = {}
if a.t1k:
    for f in glob.glob(os.path.join(a.t1k, "*", "*_genotype.tsv")):
        s = os.path.basename(f).replace("_genotype.tsv", "")
        for line in open(f):
            t = line.rstrip("\n").split("\t")
            if t[0] not in GENES or int(t[1]) == 0 or int(t[4]) <= 0: continue   # quality <= 0: spurious low-abundance call
            second = int(t[1]) > 1 and t[5] != "." and int(t[7]) > 0
            v = [norm(t[2].split(",")[0])] + ([norm(t[5].split(",")[0])] if second else [norm(t[2].split(",")[0])])
            t1k[(s, t[0])] = v

def amb_match(pub, q, k):
    """published pair with '/'-separated ambiguity lists vs a called pair, at k fields"""
    if len(q) != 2 or None in q: return np.nan
    P = [set(trunc(norm(g + "*" + x), k) for x in al.split("/")) for al, g in pub]
    Q = [trunc(x, k) for x in q]
    return float((Q[0] in P[0] and Q[1] in P[1]) or (Q[0] in P[1] and Q[1] in P[0]))
pub = {}
if a.published:
    P = pd.read_csv(a.published, sep=" ", quotechar='"', dtype=str)
    for _, r in P.iterrows():
        for g in ("A", "B", "C", "DRB1", "DQB1"):
            if isinstance(r[g], str) and isinstance(r[g + ".1"], str):
                pub[(r["id"], "HLA-" + g)] = [(r[g], "HLA-" + g), (r[g + ".1"], "HLA-" + g)]
if pub and t1k:
    prow = []
    for (s, g), pp in pub.items():
        t = t1k.get((s, g))
        if t is None: continue
        prow.append({"sample": s, "gene": g, "published": "/".join(x for x, _ in pp), "t1k": "/".join(map(str, t)),
                     "pub_t1k_1f": amb_match(pp, t, 1), "pub_t1k_2f": amb_match(pp, t, 2)})
    PR = pd.DataFrame(prow)
    if len(PR):
        PR.to_csv(a.out_prefix + "_t1k_vs_published.tsv", sep="\t", index=False)
        print("T1K vs published 1000G HLA genotypes:\n", PR.groupby("gene")[["pub_t1k_1f", "pub_t1k_2f"]].agg(["mean", "count"]).round(3).to_string())

def consensus_support(sample, gene):
    """edit distance of each FuFiHLA consensus to the closest assembled copy of the gene in the same individual
    (minimap2 asm5, NM of the best hit); FuFiHLA can split a homozygous gene into two consensuses and template one onto
    a near-identical allele, which then does not occur in the assembly"""
    cons = sorted(glob.glob(os.path.join(os.path.dirname(a.fufihla.rstrip("/")), "consensus", sample, f"{gene}_asm*.fa")))
    gf = os.path.join(a.gene_fasta_dir, f"{gene}.fa")
    if not cons or not os.path.exists(gf): return None
    with tempfile.TemporaryDirectory() as td:
        tgt = os.path.join(td, "t.fa")
        with open(tgt, "w") as out:
            keep = False
            for line in open(gf):
                if line.startswith(">"): keep = line[1:].split("#")[0] == sample
                if keep: out.write(line)
        if os.path.getsize(tgt) == 0: return None
        best = {}
        res = subprocess.run(["minimap2", "-cx", "asm5", "--secondary=no", tgt] + cons, capture_output=True, text=True).stdout
        for line in res.splitlines():
            t = line.split("\t"); nm = next(int(x[5:]) for x in t[12:] if x.startswith("NM:i:"))
            unal = int(t[1]) - (int(t[3]) - int(t[2]))
            best[t[0]] = min(best.get(t[0], 10**9), nm + unal)
        return max(best.values()) if len(best) == 2 else None

# HPRC 4-field truth (Lai et al. 2024, CSBJ): phased-assembly labels cross-checked with targeted capture sequencing
truth = {}
if a.hprc_truth:
    X = pd.read_excel(a.hprc_truth, header=None)
    hdr = X.iloc[1].tolist()
    for i in range(2, len(X)):
        smp = X.iloc[i, 1]
        if not isinstance(smp, str): continue
        for j, h in enumerate(hdr):
            if isinstance(h, str) and h.startswith("HLA-") and h in GENES:
                v = [X.iloc[i, j], X.iloc[i, j + 1]]
                v = [norm(x) for x in v if isinstance(x, str)]
                if len(v) == 1: v = v * 2
                if len(v) == 2: truth[(smp, h)] = v
    trow = []
    for (smp, g), tv in truth.items():
        r = {"sample": smp, "gene": g, "truth": "/".join(tv)}
        for meth, src in (("immuannot", imm), ("fufihla", fufi), ("t1k", t1k)):
            q = src.get((smp, g)) or src.get((smp.lower(), g))
            if q is None: continue
            r[meth] = "/".join(map(str, q))
            for k in (1, 2, 3, 4): r[f"truth_{meth}_{k}f"] = pair_match(tv, q, k)
        trow.append(r)
    TR = pd.DataFrame(trow)
    if len(TR):
        TR.to_csv(a.out_prefix + "_vs_hprc_truth.tsv", sep="\t", index=False)
        tc = [c for c in TR.columns if c.startswith("truth_")]
        print("vs HPRC 4-field truth (Lai 2024):\n", TR.groupby("gene")[tc].mean().round(3).to_string())

rows = []
for (s, g), p in imm.items():
    r = {"sample": s, "cohort": cohort.get(s), "gene": g, "immuannot": "/".join(map(str, p)), "fufihla": None, "t1k": None}
    base = s.split(".")[0]          # NA18940.JaSaPaGe -> NA18940 (1000G read data are per individual)
    q = fufi.get((s, g)); t = t1k.get((s, g)) or t1k.get((base, g)) or t1k.get((base.upper(), g))
    if q is not None:
        r["fufihla"] = "/".join(map(str, q))
        for k in (1, 2, 3, 4): r[f"imm_fufi_{k}f"] = pair_match(p, q, k)
        r["immuannot_novel"] = imm_novel.get((s, g), False); r["fufihla_inexact"] = fufi_inexact.get((s, g), False)
        r["fufihla_consensus_max_edit_to_assembly"] = consensus_support(s, g) if (a.gene_fasta_dir and r["imm_fufi_2f"] != 1.0) else None
        if r["imm_fufi_2f"] == 1.0: r["discordance_class"] = "concordant (2-field)"
        elif r["fufihla_consensus_max_edit_to_assembly"] not in (None, 0):
            # reconstructed consensus absent from the individual's assembled copies: FuFiHLA split identical (homozygous)
            # reads over two near-identical templates; genuine novel consensuses match the assembly exactly
            r["discordance_class"] = "FuFiHLA artefact (consensus not present in the assembly)"
        elif r["immuannot_novel"] and r["fufihla_inexact"]: r["discordance_class"] = "both novel (no exact IPD-IMGT/HLA match)"
        elif r["immuannot_novel"]: r["discordance_class"] = "FuFiHLA exact allele, Immuannot novel (newer IPD-IMGT/HLA release)"
        else: r["discordance_class"] = "discordant"
    if t is not None:
        r["t1k"] = "/".join(map(str, t))
        for k in (1, 2, 3): r[f"imm_t1k_{k}f"] = pair_match(p, t, k)
        if q is not None:
            for k in (1, 2, 3): r[f"fufi_t1k_{k}f"] = pair_match(q, t, k)
    rows.append(r)
R = pd.DataFrame(rows)
R.to_csv(a.out_prefix + ".tsv", sep="\t", index=False)
cols = [c for c in R.columns if re.match(r"(imm_fufi|imm_t1k|fufi_t1k)_\df", c)]
summ = R.groupby("gene")[cols].agg(["mean", "count"])
summ.to_csv(a.out_prefix + "_by_gene.tsv", sep="\t")
if "discordance_class" in R:
    dc = R.groupby(["gene", "discordance_class"]).size().unstack(fill_value=0)
    dc.to_csv(a.out_prefix + "_discordance_classes.tsv", sep="\t"); print(dc.to_string())
coh = R.groupby("cohort")[[c for c in cols if c in ("imm_fufi_2f", "imm_fufi_4f", "imm_t1k_2f")]].mean()
coh.to_csv(a.out_prefix + "_by_cohort.tsv", sep="\t")
print(summ.round(3).to_string()); print(coh.round(3).to_string())

fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))
x = np.arange(len(GENES))
for ax, (pref, title, ks) in zip(axes, [("imm_fufi", "FuFiHLA (assembly pseudo-reads) vs Immuannot", (2, 3, 4)),
                                          ("imm_t1k", "T1K (1000G Illumina reads) vs Immuannot", (1, 2, 3)),
                                          ("fufi_t1k", "T1K vs FuFiHLA", (1, 2, 3))]):
    w = 0.8 / len(ks)
    for j, k in enumerate(ks):
        c = f"{pref}_{k}f"
        if c not in R: continue
        m = R.groupby("gene")[c].mean().reindex(GENES); n = R.groupby("gene")[c].count().reindex(GENES)
        ax.bar(x + (j - (len(ks) - 1) / 2) * w, m.values, w, label=f"{k}-field")
        if j == 0:
            for xi, nn in zip(x, n.values): ax.text(xi, 1.02, f"n={int(nn) if nn == nn else 0}", ha="center", fontsize=7)
    ax.set_xticks(x); ax.set_xticklabels([g.replace("HLA-", "") for g in GENES]); ax.set_ylim(0, 1.1)
    ax.set_title(title, fontsize=10); ax.set_ylabel("fraction of individuals with identical genotype"); ax.get_legend_handles_labels()[0] and ax.legend(fontsize=8, loc="lower left")
fig.tight_layout(); fig.savefig(a.fig, dpi=150)
