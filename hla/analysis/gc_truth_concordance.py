#!/usr/bin/env python3
"""
Validate t1k, leechuck_t1k, and Immuannot (old/new DB) against real published
truth genotypes (Gourraud et al. 2014, 1000G HLA diversity panel -- Sanger/SSO
typing, 5 genes: A/B/C/DRB1/DQB1), to answer: does the pangenome-consensus
Immuannot arm improve typing accuracy over direct-read t1k, or not?

Ambiguity-aware matching (reused from leechuck's analysis/typing_concordance.py
amb_match): the published file gives, per haplotype slot, a '/'-separated list
of allele names the Sanger-era method couldn't distinguish -- a call matches
that slot if it equals ANY member of the list, at 2-field resolution.

Writes gc_truth_concordance.tsv and fig6_truth_concordance.png under
mhc_giraffe/compare/.
"""
import re
import csv
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).parent))
from gc_typing_concordance import (
    load_t1k, load_immuannot, field2, best_match_count,
    T1K_DIR, IMMU_OLD_DIR, IMMU_NEW_DIR, LEECHUCK_T1K_DIR,
    COLOR_T1K, COLOR_IMMU_OLD, COLOR_IMMU_NEW, COLOR_LEECHUCK_T1K,
    INK, INK_SECONDARY, GRID, SURFACE, style_ax,
)

W = Path("/lustre10/home/dawnxchen/mhc_giraffe")
OUT = W / "compare"
TRUTH_FILE = Path("/lustre10/home/dawnxchen/pangenome-bh26/1000g_ground_truth/1000G_2014/20140702_hla_diversity.txt")
SAMPLES = [l.strip() for l in open(W / "samples.txt") if l.strip()]
TRUTH_GENES = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1"]


def parse_truth_row(row):
    """One header line's tokens, honoring quoted fields, from a space-delimited
    quoted file (R write.table default): "id" "sbgroup" "A" "A.1" ..."""
    return re.findall(r'"([^"]*)"', row)


def load_truth():
    """Return {(sample, gene): (set_of_2field_candidates_slot1, set_..._slot2)}."""
    with open(TRUTH_FILE) as fh:
        header = parse_truth_row(fh.readline())
        out = {}
        for line in fh:
            vals = parse_truth_row(line)
            row = dict(zip(header, vals))
            s = row["id"]
            for g in ("A", "B", "C", "DRB1", "DQB1"):
                a, b = row.get(g), row.get(g + ".1")
                if not a or not b:
                    continue
                gene = "HLA-" + g
                slot1 = {field2(f"{gene}*{x}") for x in a.split("/")}
                slot2 = {field2(f"{gene}*{x}") for x in b.split("/")}
                out[(s, gene)] = (slot1 - {None}, slot2 - {None})
    return out


def truth_match(call, truth):
    """call = (a1, a2) raw allele strings (already field2-truncatable);
    truth = (set1, set2) of acceptable 2-field names per slot. Returns
    2 (both alleles matched, in either pairing), 1, or 0."""
    a1, a2 = field2(call[0]), field2(call[1])
    if a1 is None or a2 is None:
        return None
    set1, set2 = truth
    straight = (a1 in set1) + (a2 in set2)
    swapped = (a1 in set2) + (a2 in set1)
    return max(straight, swapped)


truth = load_truth()
print(f"published truth: {len(truth)} sample-gene entries")

t1k_calls, leechuck_calls, immu_old_calls, immu_new_calls = {}, {}, {}, {}
for s in SAMPLES:
    t1k = load_t1k(s, T1K_DIR)
    if t1k is not None:
        t1k_calls[s] = t1k
    lk = load_t1k(s, LEECHUCK_T1K_DIR)
    if lk is not None:
        leechuck_calls[s] = lk
    old = load_immuannot(s, IMMU_OLD_DIR, per_sample_subdir=False)
    if old is not None:
        immu_old_calls[s] = old
    new = load_immuannot(s, IMMU_NEW_DIR, per_sample_subdir=True)
    if new is not None:
        immu_new_calls[s] = new

DATASETS = [
    ("t1k", t1k_calls, COLOR_T1K),
    ("leechuck_t1k", leechuck_calls, COLOR_LEECHUCK_T1K),
    ("immu_old", immu_old_calls, COLOR_IMMU_OLD),
    ("immu_new", immu_new_calls, COLOR_IMMU_NEW),
]

# per-dataset, per-gene tallies
rows = []
pooled = {}
for name, calls, _ in DATASETS:
    n = geno_match = allele_match_sum = 0
    per_gene = defaultdict(lambda: [0, 0, 0])  # n, geno_match, allele_match_sum
    for s, gene_dict in calls.items():
        for gene in TRUTH_GENES:
            key = (s, gene)
            if key not in truth or gene not in gene_dict:
                continue
            m = truth_match(gene_dict[gene], truth[key])
            if m is None:
                continue
            n += 1
            geno_match += int(m == 2)
            allele_match_sum += m
            rec = per_gene[gene]
            rec[0] += 1
            rec[1] += int(m == 2)
            rec[2] += m
    pooled[name] = {"n": n, "geno_pct": 100 * geno_match / n if n else 0,
                     "allele_pct": 100 * allele_match_sum / (2 * n) if n else 0}
    for gene, (gn, ggm, gam) in per_gene.items():
        rows.append([name, gene, gn, f"{100*gam/(2*gn):.1f}" if gn else "0",
                     f"{100*ggm/gn:.1f}" if gn else "0"])

with open(OUT / "truth_concordance.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["dataset", "gene", "n_samples", "allele_concordance_pct", "genotype_concordance_pct"])
    for r in sorted(rows):
        w.writerow(r)
print("wrote:", OUT / "truth_concordance.tsv")

with open(OUT / "truth_concordance_overall.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["dataset", "n_gene_sample_obs", "allele_concordance_pct", "genotype_concordance_pct"])
    for name, _, _ in DATASETS:
        rec = pooled[name]
        w.writerow([name, rec["n"], f"{rec['allele_pct']:.1f}", f"{rec['geno_pct']:.1f}"])
        print(f"  {name}: n={rec['n']}  allele={rec['allele_pct']:.1f}%  genotype={rec['geno_pct']:.1f}%")
print("wrote:", OUT / "truth_concordance_overall.tsv")

# --- Figure: concordance vs published truth, pooled over 5 genes ---
fig, ax = plt.subplots(figsize=(6.5, 4.5), facecolor=SURFACE)
labels = ["t1k\n(your run)", "leechuck_t1k\n(independent)", "Immuannot\n(old DB)", "Immuannot\n(new DB)"]
values = [pooled[name]["geno_pct"] for name, _, _ in DATASETS]
colors = [c for _, _, c in DATASETS]
ns = [pooled[name]["n"] for name, _, _ in DATASETS]
bars = ax.bar(labels, values, color=colors, width=0.6, zorder=3)
for b, v, n in zip(bars, values, ns):
    ax.text(b.get_x() + b.get_width()/2, v + 1.5, f"{v:.0f}%\n(n={n})", ha="center", va="bottom", fontsize=9, color=INK)
ax.set_ylim(0, 105)
ax.set_ylabel("full-genotype concordance (%) vs. published truth\n(Gourraud 2014, 2-field, 5 genes)", color=INK_SECONDARY, fontsize=9)
ax.set_title("Which typing result matches real published genotypes?", color=INK, fontsize=12, loc="left", pad=12)
style_ax(ax)
fig.tight_layout()
fig.savefig(OUT / "fig6_truth_concordance.png", dpi=150, facecolor=SURFACE)
plt.close(fig)
print("wrote:", OUT / "fig6_truth_concordance.png")
