#!/usr/bin/env python3
"""Figure 5: HLA allele landscape by cohort (allele-frequency heatmaps, DR haplogroups, C4 forms, novel-allele fraction).
Run inside hla-viz (hla_calls.tsv, gene_copy_number.tsv from the aggregate step)."""
import pandas as pd, numpy as np, matplotlib
from cohorts import ORDER, ONE, COLOR, present
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d=pd.read_csv("hla_calls.tsv",sep="\t",low_memory=False); d=d[d.cohort!="REF"]
cn=pd.read_csv("gene_copy_number.tsv",sep="\t"); cn=cn[cn.cohort!="REF"]
order=present(d.cohort)
n={c:d[d.cohort==c].hap_id.nunique() for c in order}
short={c:ONE[c].replace(" (","\n(")+f"\nn={n[c]}" for c in order}
one=ONE
fig=plt.figure(figsize=(22,12.5)); gs=fig.add_gridspec(2,3,height_ratios=[3,1.3],hspace=0.5,wspace=0.35)
axes=[fig.add_subplot(gs[0,i]) for i in range(3)]+[fig.add_subplot(gs[1,i]) for i in range(3)]
def freq_table(gene):
    s=d[d.gene==gene]
    f=s.groupby("cohort").allele_2field.value_counts(normalize=True).unstack(0).reindex(columns=order).fillna(0)
    keep=set()
    for c in order: keep|=set(f[c].sort_values(ascending=False).head(6).index)
    f=f.loc[sorted(keep,key=lambda a:-f.loc[a].max())]; f.index=[a.replace(gene+"*","") for a in f.index]; return f
for ax,gene in zip(axes[:3],["HLA-A","HLA-B","HLA-DRB1"]):
    f=freq_table(gene); im=ax.imshow(f.values,cmap="YlOrRd",vmin=0,vmax=0.35,aspect="auto")
    ax.set_yticks(range(len(f))); ax.set_yticklabels(f.index,fontsize=9); ax.set_xticks(range(len(order))); ax.set_xticklabels([short[c] for c in order],fontsize=8); ax.set_title(f"{gene} two-field alleles",fontsize=11)
    for i in range(len(f)):
        for j in range(len(order)):
            v=f.values[i,j]
            if v>0: ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=6,color="white" if v>0.2 else "black")
fig.colorbar(im,ax=axes[:3],fraction=0.015,pad=0.01,label="haplotype frequency")
ax=axes[3]; drb=cn.groupby("cohort")[["HLA-DRB3","HLA-DRB4","HLA-DRB5"]].mean().reindex(order); drb["none"]=1-drb.sum(axis=1)
drb.clip(lower=0).plot.bar(stacked=True,ax=ax,color=["#1f77b4","#ff7f0e","#2ca02c","#cccccc"],width=0.7)
ax.set_xticklabels([one[c] for c in order],rotation=45,ha="right",fontsize=9); ax.set_ylabel("fraction of haplotypes"); ax.set_title("Secondary DRB gene (DR haplogroup)",fontsize=10); ax.legend(fontsize=7,ncol=2,loc="lower left"); ax.set_xlabel("")
ax=axes[4]; c4=cn.groupby("cohort")[["C4AL","C4AS","C4BL","C4BS"]].mean().reindex(order)
c4.plot.bar(ax=ax,width=0.8,color=["#d62728","#ff9896","#9467bd","#c5b0d5"]); ax.set_xticklabels([one[c] for c in order],rotation=45,ha="right",fontsize=9); ax.set_ylabel("mean copies per haplotype"); ax.set_title("C4A / C4B, long (L) and short (S) forms",fontsize=10); ax.legend(fontsize=7,ncol=2); ax.set_xlabel("")
ax=axes[5]; cl=["HLA-A","HLA-B","HLA-C","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1"]
x=d[d.gene.isin(cl)].copy(); x["novel"]=x.novel.astype(str)=="True"
nv=x.groupby(["gene","cohort"]).novel.mean().unstack("cohort").reindex(index=cl,columns=order)
nv.plot.bar(ax=ax,width=0.85,color=[COLOR[c] for c in order]); ax.set_xticklabels([g.replace("HLA-","") for g in cl],rotation=0,fontsize=8); ax.set_ylabel("fraction of gene copies"); ax.set_title("Full-length allele absent from IPD-IMGT/HLA 3.55",fontsize=10); ax.legend([one[c] for c in order],fontsize=6.5,ncol=3); ax.set_xlabel("")
fig.suptitle(f"HLA allele landscape by cohort (HPRC r2 split by population; HPRC Jewish = HG002 only) across {sum(n.values())} assembled haplotypes (Immuannot on assembly-derived MHC)",fontsize=13,y=0.98)
fig.savefig("figures/fig5_population_hla.png",dpi=160,bbox_inches="tight")
