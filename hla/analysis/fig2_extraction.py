#!/usr/bin/env python3
"""Figure 2: MHC extraction completeness per cohort (fraction of the GRCh38 MHC covered by the extracted segments,
summed over segments of a haplotype) and number of segments. Run inside the viz dir (mhc_extraction_summary.tsv)."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from cohorts import ONE, COLOR, present
m = pd.read_csv("mhc_extraction_summary.tsv", sep="\t"); m = m[m.cohort != "REF"]
h = m.groupby(["cohort", "sample", "haplotype"]).agg(cov=("mhc_covered_frac", "sum"), pieces=("n_pieces", "first")).reset_index()
h["cov"] = h["cov"].clip(upper=1.0)
order = present(h.cohort)
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [2.2, 1]})
rng = np.random.default_rng(0)
for i, c in enumerate(order):
    v = h[h.cohort == c]["cov"].values
    ax.scatter(i + rng.uniform(-0.25, 0.25, len(v)), v, s=9, color=COLOR[c], alpha=0.7, lw=0)
    ax.text(i, 0.905, f"{(v >= 0.99).sum()}/{len(v)}", ha="center", fontsize=8)
ax.set_xticks(range(len(order))); ax.set_xticklabels([ONE[c] for c in order], rotation=35, ha="right", fontsize=9)
ax.set_ylim(0.9, 1.01); ax.set_ylabel("fraction of GRCh38 MHC covered")
ax.set_title("MHC coverage per haplotype (numbers: haplotypes with coverage >= 0.99)", fontsize=10, loc="left")
pc = pd.crosstab(h.cohort, h.pieces.clip(upper=3)).reindex(order)
pc = pc.div(pc.sum(axis=1), axis=0)
bottom = np.zeros(len(order))
for k, col in zip(pc.columns, ["#4c72b0", "#dd8452", "#c44e52"]):
    ax2.bar(range(len(order)), pc[k].values, bottom=bottom, color=col, label={1: "1 segment", 2: "2 segments", 3: "3+ segments"}[k]); bottom += pc[k].values
ax2.set_xticks(range(len(order))); ax2.set_xticklabels([ONE[c] for c in order], rotation=35, ha="right", fontsize=9)
ax2.set_ylabel("fraction of haplotypes"); ax2.legend(fontsize=8, loc="lower right"); ax2.set_title("Segments per haplotype", fontsize=10, loc="left")
fig.tight_layout(); fig.savefig("figures/fig2_mhc_extraction.png", dpi=170)
print("haplotypes", len(h), "coverage>=0.99", int((h["cov"] >= 0.99).sum()), "single segment", int((h.pieces == 1).sum()))
print(h.groupby("cohort").agg(n=("cov", "size"), cov99=("cov", lambda v: int((v >= 0.99).sum())), single=("pieces", lambda v: int((v == 1).sum()))).reindex(order).to_string())
