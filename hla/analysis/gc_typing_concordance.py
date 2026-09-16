#!/usr/bin/env python3
"""
Compare HLA/MHC typing calls across four pipelines run on (subsets of) the same
1000G subsampled short-read cohort:

  t1k          - run-t1k directly on raw FASTQs (mhc_giraffe/t1k/<S>/<S>_genotype.tsv)
  immu_old     - Immuannot on phased pangenome consensus, IPD-IMGT/HLA 3.55.0 db
                 (mhc_giraffe/immuannot/<S>.hap{1,2}.gtf.gz)
  immu_new     - Immuannot on the same consensus, IPD-IMGT/HLA 3.65.0 db
                 (mhc_giraffe/immuannot_newdb/<S>/<S>.hap{1,2}.gtf.gz)
  leechuck_t1k - leechuck's own run-t1k on raw 1000G Illumina reads, no pangenome
                 involved at all (leechuck/hla/typing/conc_partial/data/
                 conc_t1k_vs_published.tsv, 'sample'/'gene'/'t1k' columns).
                 Only covers 5 classical genes (A/B/C/DQB1/DRB1) x 412 samples
                 (whichever samples had a published truth genotype to merge
                 against) -- NOT the same denominator as the other 3 datasets,
                 so it is reported separately (fig5) rather than folded into
                 the whole-cohort completion/concordance figures.

Writes summary TSVs and PNG figures under mhc_giraffe/compare/.
"""
import gzip
import re
import csv
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = Path("/lustre10/home/dawnxchen/mhc_giraffe")
T1K_DIR = W / "t1k"
IMMU_OLD_DIR = W / "immuannot"
IMMU_NEW_DIR = W / "immuannot_newdb"
LEECHUCK_T1K_TSV = Path("/lustre10/home/leechuck/hla/typing/conc_partial/data/conc_t1k_vs_published.tsv")
OUT = W / "compare"
OUT.mkdir(exist_ok=True)

SAMPLES = [l.strip() for l in open(W / "samples.txt") if l.strip()]

# Classical HLA genes typed by both t1k and Immuannot; focus of the per-gene plots.
CLASSICAL_GENES = [
    "HLA-A", "HLA-B", "HLA-C",
    "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
]
LEECHUCK_T1K_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DQB1", "HLA-DRB1"]

# validated categorical palette (dataviz skill, references/palette.md), slots 1-4
COLOR_T1K = "#2a78d6"          # blue
COLOR_IMMU_OLD = "#eb6834"     # orange
COLOR_IMMU_NEW = "#1baf7a"     # aqua
COLOR_LEECHUCK_T1K = "#eda100"  # yellow
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e1e0d9"
SURFACE = "#fcfcfb"


def field1(allele: str):
    """Truncate to 1-field (allele-group) resolution, e.g.
    'HLA-A*29:02:01:01' -> 'HLA-A*29', 'HLA-A*01:new' -> 'HLA-A*01'."""
    if not allele or allele == ".":
        return None
    a = allele[:-4] if allele.endswith(":new") else allele
    return a.split(":")[0]


def field2(allele: str):
    """Truncate to 2-field resolution, e.g. 'HLA-A*29:02:01:01' -> 'HLA-A*29:02'.
    Immuannot's ':new' suffix marks the point past which it couldn't resolve a
    known allele -- e.g. 'HLA-A*01:new' is only 1-field-resolved, NOT a literal
    2nd field of 'new'. Strip it first, then only return a 2-field truncation if
    at least 2 real fields remain; otherwise this allele isn't comparable at
    2-field resolution."""
    if not allele or allele == ".":
        return None
    a = allele[:-4] if allele.endswith(":new") else allele
    parts = a.split(":")
    if len(parts) < 2:
        return None
    return ":".join(parts[:2])


def load_t1k(sample: str):
    """Return {gene: (allele1, allele2)} raw strings, or None if not run yet."""
    f = T1K_DIR / sample / f"{sample}_genotype.tsv"
    if not f.exists() or f.stat().st_size == 0:
        return None
    out = {}
    with open(f) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            gene, _n, a1, _s1, _ab1, a2 = cols[:6]
            out[gene] = (a1 if a1 != "." else None, a2 if a2 != "." else None)
    return out


GENE_RE = re.compile(r'gene_name "([^"]+)"')
CONSENSUS_RE = re.compile(r'consensus "([^"]+)"')


def load_immuannot_hap(gtf_gz: Path):
    """Return {gene: consensus_allele} from one hap*.gtf.gz transcript line."""
    out = {}
    if not gtf_gz.exists() or gtf_gz.stat().st_size == 0:
        return out
    with gzip.open(gtf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            gm = GENE_RE.search(cols[8])
            cm = CONSENSUS_RE.search(cols[8])
            if gm and cm:
                out[gm.group(1)] = cm.group(1)
    return out


def load_immuannot(sample: str, base: Path, per_sample_subdir: bool):
    """Return {gene: (hap1_allele, hap2_allele)} or None if not run yet."""
    d = base / sample if per_sample_subdir else base
    f1 = d / f"{sample}.hap1.gtf.gz"
    f2 = d / f"{sample}.hap2.gtf.gz"
    if not f1.exists() or not f2.exists():
        return None
    h1 = load_immuannot_hap(f1)
    h2 = load_immuannot_hap(f2)
    genes = set(h1) | set(h2)
    return {g: (h1.get(g), h2.get(g)) for g in genes}


def is_novel(allele: str) -> bool:
    return bool(allele) and allele.endswith(":new")


def load_leechuck_t1k():
    """Return {sample: {gene: (allele1, allele2)}} from leechuck's own
    t1k-on-raw-1000G-reads column in conc_t1k_vs_published.tsv."""
    out = defaultdict(dict)
    with open(LEECHUCK_T1K_TSV) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            alleles = row["t1k"].split("/")
            a1 = alleles[0] if alleles[0] not in ("", ".") else None
            a2 = alleles[1] if len(alleles) > 1 and alleles[1] not in ("", ".") else None
            out[row["sample"]][row["gene"]] = (a1, a2)
    return dict(out)


# ---------------------------------------------------------------------------
# Pass 1: load everything that's completed so far, track completion counts
# ---------------------------------------------------------------------------
t1k_calls, immu_old_calls, immu_new_calls = {}, {}, {}

for s in SAMPLES:
    t1k = load_t1k(s)
    if t1k is not None:
        t1k_calls[s] = t1k
    old = load_immuannot(s, IMMU_OLD_DIR, per_sample_subdir=False)
    if old is not None:
        immu_old_calls[s] = old
    new = load_immuannot(s, IMMU_NEW_DIR, per_sample_subdir=True)
    if new is not None:
        immu_new_calls[s] = new

leechuck_t1k_calls = load_leechuck_t1k()

n_total = len(SAMPLES)
completion = {
    "t1k": len(t1k_calls),
    "immu_old": len(immu_old_calls),
    "immu_new": len(immu_new_calls),
}
print(f"total samples in cohort: {n_total}")
for k, v in completion.items():
    print(f"  {k}: {v} completed ({100*v/n_total:.1f}%)")
print(f"  leechuck_t1k: {len(leechuck_t1k_calls)} samples, "
      f"{len(LEECHUCK_T1K_GENES)} genes only (separate scope, see fig5)")

# ---------------------------------------------------------------------------
# Pass 2: pairwise concordance per gene, at both 1-field (allele-group) and
# 2-field (protein) resolution.
# gene -> pair -> {"n": samples with both typed, "allele_match": sum of 0-2,
#                  "geno_match": full 2/2 matches}
# ---------------------------------------------------------------------------

RESOLUTIONS = {"field1": field1, "field2": field2}


def geno_at(gene_dict, gene, trunc_fn):
    if gene_dict is None or gene not in gene_dict:
        return None
    a, b = gene_dict[gene]
    a_t, b_t = trunc_fn(a), trunc_fn(b)
    if a_t is None and b_t is None:
        return None
    return (a_t, b_t)


def best_match_count(g1, g2):
    """Count of alleles in common between two unordered genotype pairs,
    trying both phasings, ignoring None slots."""
    a1, a2 = g1
    b1, b2 = g2
    straight = (a1 is not None and a1 == b1) + (a2 is not None and a2 == b2)
    swapped = (a1 is not None and a1 == b2) + (a2 is not None and a2 == b1)
    return max(straight, swapped)


PAIRS = [
    ("t1k", "immu_old", t1k_calls, immu_old_calls),
    ("t1k", "immu_new", t1k_calls, immu_new_calls),
    ("immu_old", "immu_new", immu_old_calls, immu_new_calls),
]

all_genes = set()
for d in (t1k_calls, immu_old_calls, immu_new_calls):
    for s, g in d.items():
        all_genes.update(g.keys())

# concordance[resolution][(name_a, name_b)][gene] = {...}
concordance = {
    res: defaultdict(lambda: defaultdict(lambda: {"n": 0, "allele_match": 0, "geno_match": 0}))
    for res in RESOLUTIONS
}

for name_a, name_b, calls_a, calls_b in PAIRS:
    common_samples = set(calls_a) & set(calls_b)
    for res, trunc_fn in RESOLUTIONS.items():
        for gene in sorted(all_genes):
            for s in common_samples:
                ga = geno_at(calls_a[s], gene, trunc_fn)
                gb = geno_at(calls_b[s], gene, trunc_fn)
                if ga is None or gb is None:
                    continue
                m = best_match_count(ga, gb)
                rec = concordance[res][(name_a, name_b)][gene]
                rec["n"] += 1
                rec["allele_match"] += m
                rec["geno_match"] += int(m == 2)

with open(OUT / "concordance_per_gene.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["resolution", "pair", "gene", "n_samples_both_typed", "allele_concordance_pct", "genotype_concordance_pct"])
    for res in RESOLUTIONS:
        for (name_a, name_b), genes in concordance[res].items():
            for gene, rec in sorted(genes.items()):
                if rec["n"] == 0:
                    continue
                ac = 100 * rec["allele_match"] / (2 * rec["n"])
                gc = 100 * rec["geno_match"] / rec["n"]
                w.writerow([res, f"{name_a}_vs_{name_b}", gene, rec["n"], f"{ac:.1f}", f"{gc:.1f}"])

# overall (all genes pooled) concordance per pair, per resolution
overall = {res: {} for res in RESOLUTIONS}
for res in RESOLUTIONS:
    for (name_a, name_b), genes in concordance[res].items():
        n = sum(r["n"] for r in genes.values())
        am = sum(r["allele_match"] for r in genes.values())
        gm = sum(r["geno_match"] for r in genes.values())
        overall[res][(name_a, name_b)] = {
            "n": n,
            "allele_pct": 100 * am / (2 * n) if n else 0,
            "geno_pct": 100 * gm / n if n else 0,
        }

with open(OUT / "concordance_overall.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["resolution", "pair", "n_gene_sample_obs", "allele_concordance_pct", "genotype_concordance_pct"])
    for res in RESOLUTIONS:
        for (a, b), rec in overall[res].items():
            w.writerow([res, f"{a}_vs_{b}", rec["n"], f"{rec['allele_pct']:.1f}", f"{rec['geno_pct']:.1f}"])

# ---------------------------------------------------------------------------
# Pass 3: novel/uncalled-allele rate per dataset (classical genes)
# ---------------------------------------------------------------------------
novel_rate = {"t1k_untyped": defaultdict(lambda: [0, 0]),
              "immu_old_new": defaultdict(lambda: [0, 0]),
              "immu_new_new": defaultdict(lambda: [0, 0])}

for s, g in t1k_calls.items():
    for gene in CLASSICAL_GENES:
        if gene not in g:
            continue
        a1, a2 = g[gene]
        cnt = novel_rate["t1k_untyped"][gene]
        cnt[1] += 2
        cnt[0] += (a1 is None) + (a2 is None)

for s, g in immu_old_calls.items():
    for gene in CLASSICAL_GENES:
        if gene not in g:
            continue
        a1, a2 = g[gene]
        cnt = novel_rate["immu_old_new"][gene]
        cnt[1] += 2
        cnt[0] += is_novel(a1) + is_novel(a2)

for s, g in immu_new_calls.items():
    for gene in CLASSICAL_GENES:
        if gene not in g:
            continue
        a1, a2 = g[gene]
        cnt = novel_rate["immu_new_new"][gene]
        cnt[1] += 2
        cnt[0] += is_novel(a1) + is_novel(a2)

with open(OUT / "novel_uncalled_rate.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["dataset", "gene", "flagged", "total_alleles", "pct"])
    for dsname, genes in novel_rate.items():
        for gene, (flagged, total) in sorted(genes.items()):
            pct = 100 * flagged / total if total else 0
            w.writerow([dsname, gene, flagged, total, f"{pct:.1f}"])

print("wrote:", OUT / "concordance_per_gene.tsv")
print("wrote:", OUT / "concordance_overall.tsv")
print("wrote:", OUT / "novel_uncalled_rate.tsv")

# ---------------------------------------------------------------------------
# Pass 4: leechuck_t1k vs the other 3 datasets, restricted to the 5 genes and
# 412 samples leechuck_t1k actually covers (different scope from Pass 2/3,
# see module docstring) -- field2 (2-field / protein) resolution only.
# ---------------------------------------------------------------------------
LEECHUCK_PAIRS = [
    ("leechuck_t1k", "t1k", leechuck_t1k_calls, t1k_calls),
    ("leechuck_t1k", "immu_old", leechuck_t1k_calls, immu_old_calls),
    ("leechuck_t1k", "immu_new", leechuck_t1k_calls, immu_new_calls),
]

leechuck_concordance = defaultdict(lambda: defaultdict(lambda: {"n": 0, "allele_match": 0, "geno_match": 0}))
for name_a, name_b, calls_a, calls_b in LEECHUCK_PAIRS:
    common_samples = set(calls_a) & set(calls_b)
    for gene in LEECHUCK_T1K_GENES:
        for s in common_samples:
            ga = geno_at(calls_a[s], gene, field2)
            gb = geno_at(calls_b[s], gene, field2)
            if ga is None or gb is None:
                continue
            m = best_match_count(ga, gb)
            rec = leechuck_concordance[(name_a, name_b)][gene]
            rec["n"] += 1
            rec["allele_match"] += m
            rec["geno_match"] += int(m == 2)

with open(OUT / "concordance_leechuck_t1k.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["pair", "gene", "n_samples_both_typed", "allele_concordance_pct", "genotype_concordance_pct"])
    for (name_a, name_b), genes in leechuck_concordance.items():
        for gene, rec in sorted(genes.items()):
            if rec["n"] == 0:
                continue
            ac = 100 * rec["allele_match"] / (2 * rec["n"])
            gc = 100 * rec["geno_match"] / rec["n"]
            w.writerow([f"{name_a}_vs_{name_b}", gene, rec["n"], f"{ac:.1f}", f"{gc:.1f}"])
print("wrote:", OUT / "concordance_leechuck_t1k.tsv")

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def style_ax(ax):
    ax.set_facecolor(SURFACE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.spines["bottom"].set_color("#c3c2b7")
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)


# --- Figure 1: completion status ---
fig, ax = plt.subplots(figsize=(6, 4), facecolor=SURFACE)
labels = ["t1k\n(raw reads)", "Immuannot\n(old DB)", "Immuannot\n(new DB)"]
values = [completion["t1k"], completion["immu_old"], completion["immu_new"]]
colors = [COLOR_T1K, COLOR_IMMU_OLD, COLOR_IMMU_NEW]
bars = ax.bar(labels, values, color=colors, width=0.6, zorder=3)
ax.axhline(n_total, color=INK_SECONDARY, linewidth=1, linestyle="--", zorder=2)
ax.text(2.55, n_total, f" cohort = {n_total}", va="center", ha="left",
        color=INK_SECONDARY, fontsize=8)
for b, v in zip(bars, values):
    ax.text(b.get_x() + b.get_width() / 2, v + n_total * 0.015, f"{v}\n({100*v/n_total:.0f}%)",
            ha="center", va="bottom", fontsize=9, color=INK)
ax.set_ylim(0, n_total * 1.15)
ax.set_ylabel("samples completed", color=INK_SECONDARY, fontsize=9)
ax.set_title("Pipeline completion so far", color=INK, fontsize=12, loc="left", pad=12)
style_ax(ax)
fig.tight_layout()
fig.savefig(OUT / "fig1_completion.png", dpi=150, facecolor=SURFACE)
plt.close(fig)

# --- Figure 2: overall full-genotype concordance, field1 vs field2 resolution, by pair ---
fig, ax = plt.subplots(figsize=(6, 4), facecolor=SURFACE)
pair_labels = ["t1k vs\nImmuannot-old", "t1k vs\nImmuannot-new", "Immuannot-old vs\nImmuannot-new"]
pair_keys = [("t1k", "immu_old"), ("t1k", "immu_new"), ("immu_old", "immu_new")]
geno_pct_f1 = [overall["field1"].get(k, {"geno_pct": 0})["geno_pct"] for k in pair_keys]
geno_pct_f2 = [overall["field2"].get(k, {"geno_pct": 0})["geno_pct"] for k in pair_keys]
x = range(len(pair_labels))
w_ = 0.32
b1 = ax.bar([i - w_/2 for i in x], geno_pct_f1, width=w_, label="1-field (allele group)", color=COLOR_T1K, zorder=3)
b2 = ax.bar([i + w_/2 for i in x], geno_pct_f2, width=w_, label="2-field (protein)", color=COLOR_IMMU_OLD, zorder=3)
for bars in (b1, b2):
    for b in bars:
        h = b.get_height()
        ax.text(b.get_x() + b.get_width()/2, h + 1.5, f"{h:.0f}%", ha="center", va="bottom", fontsize=8, color=INK)
ax.set_xticks(list(x))
ax.set_xticklabels(pair_labels, fontsize=9, color=INK_SECONDARY)
ax.set_ylim(0, 105)
ax.set_ylabel("full-genotype (2/2) concordance (%), all genes pooled", color=INK_SECONDARY, fontsize=9)
ax.set_title("Overall typing concordance between pipelines", color=INK, fontsize=12, loc="left", pad=12)
ax.legend(frameon=False, loc="lower left", fontsize=8)
style_ax(ax)
fig.tight_layout()
fig.savefig(OUT / "fig2_overall_concordance.png", dpi=150, facecolor=SURFACE)
plt.close(fig)

# --- Figure 3: per-gene genotype concordance for classical HLA genes, 2-field ---
fig, ax = plt.subplots(figsize=(10, 5), facecolor=SURFACE)
gene_list = [g for g in CLASSICAL_GENES if any(g in concordance["field2"][k] for k in pair_keys)]
n_pairs = len(pair_keys)
bar_w = 0.8 / n_pairs
colors3 = [COLOR_T1K, COLOR_IMMU_OLD, COLOR_IMMU_NEW]
for i, (pk, lab, col) in enumerate(zip(pair_keys, pair_labels, colors3)):
    ys = []
    for gene in gene_list:
        rec = concordance["field2"][pk].get(gene, {"n": 0, "geno_match": 0})
        ys.append(100 * rec["geno_match"] / rec["n"] if rec["n"] else 0)
    xs = [j + (i - (n_pairs-1)/2) * bar_w for j in range(len(gene_list))]
    ax.bar(xs, ys, width=bar_w * 0.95, label=lab.replace("\n", " "), color=col, zorder=3)
ax.set_xticks(range(len(gene_list)))
ax.set_xticklabels(gene_list, rotation=30, ha="right", fontsize=9, color=INK_SECONDARY)
ax.set_ylim(0, 105)
ax.set_ylabel("full-genotype concordance (%)", color=INK_SECONDARY, fontsize=9)
ax.set_title("Per-gene genotype concordance, classical HLA genes (2-field)", color=INK, fontsize=12, loc="left", pad=12)
ax.legend(frameon=False, loc="lower right", fontsize=8, ncol=1)
style_ax(ax)
fig.tight_layout()
fig.savefig(OUT / "fig3_per_gene_concordance.png", dpi=150, facecolor=SURFACE)
plt.close(fig)

# --- Figure 4: novel/uncalled allele rate per dataset, classical genes ---
fig, ax = plt.subplots(figsize=(10, 5), facecolor=SURFACE)
ds_order = [("t1k_untyped", "t1k: untyped alleles", COLOR_T1K),
            ("immu_old_new", "Immuannot-old: novel (:new) calls", COLOR_IMMU_OLD),
            ("immu_new_new", "Immuannot-new: novel (:new) calls", COLOR_IMMU_NEW)]
n_ds = len(ds_order)
bar_w = 0.8 / n_ds
for i, (key, lab, col) in enumerate(ds_order):
    ys = []
    for gene in CLASSICAL_GENES:
        flagged, total = novel_rate[key].get(gene, [0, 0])
        ys.append(100 * flagged / total if total else 0)
    xs = [j + (i - (n_ds-1)/2) * bar_w for j in range(len(CLASSICAL_GENES))]
    ax.bar(xs, ys, width=bar_w * 0.95, label=lab, color=col, zorder=3)
ax.set_xticks(range(len(CLASSICAL_GENES)))
ax.set_xticklabels(CLASSICAL_GENES, rotation=30, ha="right", fontsize=9, color=INK_SECONDARY)
ax.set_ylabel("% of alleles flagged", color=INK_SECONDARY, fontsize=9)
ax.set_title("Untyped (t1k) / novel-allele (Immuannot) rate by gene", color=INK, fontsize=12, loc="left", pad=12)
ax.legend(frameon=False, loc="upper right", fontsize=8)
style_ax(ax)
fig.tight_layout()
fig.savefig(OUT / "fig4_novel_uncalled_rate.png", dpi=150, facecolor=SURFACE)
plt.close(fig)

# --- Figure 5: leechuck's raw-read t1k vs the 3 pipeline datasets, 2-field,
# restricted to the 5 genes / 412-sample scope it actually covers ---
fig, ax = plt.subplots(figsize=(8, 5), facecolor=SURFACE)
lk_pair_keys = [("leechuck_t1k", "t1k"), ("leechuck_t1k", "immu_old"), ("leechuck_t1k", "immu_new")]
lk_pair_labels = ["vs your t1k\n(raw reads)", "vs Immuannot\n(old DB)", "vs Immuannot\n(new DB)"]
lk_colors = [COLOR_T1K, COLOR_IMMU_OLD, COLOR_IMMU_NEW]
n_pairs = len(lk_pair_keys)
bar_w = 0.8 / n_pairs
for i, (pk, lab, col) in enumerate(zip(lk_pair_keys, lk_pair_labels, lk_colors)):
    ys = []
    for gene in LEECHUCK_T1K_GENES:
        rec = leechuck_concordance[pk].get(gene, {"n": 0, "geno_match": 0})
        ys.append(100 * rec["geno_match"] / rec["n"] if rec["n"] else 0)
    xs = [j + (i - (n_pairs - 1) / 2) * bar_w for j in range(len(LEECHUCK_T1K_GENES))]
    ax.bar(xs, ys, width=bar_w * 0.95, label=lab.replace("\n", " "), color=col, zorder=3)
ax.set_xticks(range(len(LEECHUCK_T1K_GENES)))
ax.set_xticklabels(LEECHUCK_T1K_GENES, fontsize=9, color=INK_SECONDARY)
ax.set_ylim(0, 105)
ax.set_ylabel("full-genotype concordance (%), 2-field", color=INK_SECONDARY, fontsize=9)
ax.set_title("leechuck's raw-read t1k vs. your 3 datasets (n up to 412, 5 genes)",
              color=INK, fontsize=12, loc="left", pad=12)
ax.legend(frameon=False, loc="lower right", fontsize=8)
style_ax(ax)
fig.tight_layout()
fig.savefig(OUT / "fig5_leechuck_t1k_comparison.png", dpi=150, facecolor=SURFACE)
plt.close(fig)

print("wrote figures to", OUT)
