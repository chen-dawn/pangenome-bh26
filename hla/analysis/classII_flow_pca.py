#!/usr/bin/env python3
"""Class II gene-level haplotype strings, alluvial haplotype-flow figures (all haplotypes; per cohort) and a
diplotype PCA on two-field class II alleles - after Chin, 'Human pangenome graph analysis for the MHC' (ASHI 2023).
Run inside hla-viz (hla_calls.tsv)."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.path import Path
from cohorts import ORDER, SHORT, COLOR, present
d=pd.read_csv("hla_calls.tsv",sep="\t",low_memory=False); d=d[d.cohort!="REF"]
genes=["HLA-DRA","DRB345","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DQA2","HLA-DQB2","HLA-DOB","TAP2","TAP1","HLA-DMB","HLA-DMA"]
def two(a): return "-" if pd.isna(a) else a.replace("HLA-","")
hap={}
for h,s in d.groupby("hap_id"):
    row={"cohort":s.cohort.iloc[0]}
    for g in genes:
        if g=="DRB345":
            x=s[s.gene.isin(["HLA-DRB3","HLA-DRB4","HLA-DRB5"])]; row[g]="none" if x.empty else "/".join(sorted(set(two(a) for a in x.allele_2field)))
        else:
            x=s[s.gene==g]; row[g]="-" if x.empty else two(x.allele_2field.iloc[0])
    hap[h]=row
H=pd.DataFrame(hap).T; H.to_csv("data/classII_haplotype_strings.tsv",sep="\t")
H["string"]=H[genes].agg("-".join,axis=1)
top=H.string.value_counts(); print("distinct class II gene-level haplotypes:",len(top),"of",len(H)); print(top.head(12).to_string())
def alluvial(ax,H,cols,color_by,title,min_n=1):
    xp=np.arange(len(cols))*1.0; n=len(H); order={}; ypos={}
    for j,c in enumerate(cols):
        vc=H[c].value_counts(); vc=vc[vc>=min_n]; order[c]=list(vc.index); y=0; ypos[c]={}; gap=0.006*n
        for k in order[c]: ypos[c][k]=(y,y+vc[k]); y+=vc[k]+gap
        for k in order[c]:
            y0,y1=ypos[c][k]; ax.add_patch(plt.Rectangle((xp[j]-0.06,y0),0.12,y1-y0,color="#444444",lw=0))
            if y1-y0>0.012*n: ax.text(xp[j]+0.08 if j==len(cols)-1 else xp[j]-0.08,(y0+y1)/2,k,fontsize=6.5,va="center",ha="left" if j==len(cols)-1 else "right")
    groups=sorted(H[color_by].unique()); cmap=plt.get_cmap("tab20"); colr={g:cmap(i%20) for i,g in enumerate(groups)}
    for j in range(len(cols)-1):
        a,b=cols[j],cols[j+1]; cur_a={k:ypos[a][k][0] for k in ypos[a]}; cur_b={k:ypos[b][k][0] for k in ypos[b]}
        flows=H.groupby([a,b,color_by]).size()
        for k in order[a]:
            for kb in order[b]:
                for g in groups:
                    m=flows.get((k,kb,g),0)
                    if not m: continue
                    ya=cur_a[k]; yb=cur_b[kb]
                    verts=[(xp[j]+0.06,ya),(xp[j]+0.5,ya),(xp[j+1]-0.5,yb),(xp[j+1]-0.06,yb),(xp[j+1]-0.06,yb+m),(xp[j+1]-0.5,yb+m),(xp[j]+0.5,ya+m),(xp[j]+0.06,ya+m),(xp[j]+0.06,ya)]
                    codes=[Path.MOVETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.LINETO,Path.CURVE4,Path.CURVE4,Path.CURVE4,Path.CLOSEPOLY]
                    ax.add_patch(PathPatch(Path(verts,codes),facecolor=colr[g],alpha=0.6,lw=0)); cur_a[k]+=m; cur_b[kb]+=m
    ax.set_xlim(-0.6,len(cols)-0.4); ax.set_ylim(-0.01*n,n*1.08); ax.set_xticks(xp); ax.set_xticklabels([c.replace("HLA-","") for c in cols]); ax.set_yticks([]); ax.set_title(title,fontsize=10,loc="left")
    for s in ["top","right","left"]: ax.spines[s].set_visible(False)
    return colr
H["DRB1grp"]=H["HLA-DRB1"].str.split(":").str[0]
H["DRhap"]=H["DRB345"].map(lambda v: "DRB3" if v.startswith("DRB3") else "DRB4" if v.startswith("DRB4") else "DRB5" if v.startswith("DRB5") else "none")
cols=["DRB345","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DQA2","HLA-DQB2","TAP2","TAP1"]
fig,ax=plt.subplots(figsize=(15,9))
colr=alluvial(ax,H,cols,"DRhap",f"Gene-level haplotype flow, {len(H)} haplotypes (ribbons coloured by secondary DRB gene: DRB3 / DRB4 / DRB5 / none)")
ax.legend(handles=[plt.Rectangle((0,0),1,1,color=colr[g],alpha=0.6) for g in colr],labels=list(colr),fontsize=8,loc="upper right")
fig.tight_layout(); fig.savefig("figures/fig6_classII_haplotype_flow.png",dpi=160)
order=present(H.cohort)
fig,axes=plt.subplots(2,(len(order)+1)//2,figsize=(6*((len(order)+1)//2),16)); axes=axes.ravel()
for ax in axes[len(order):]: ax.axis("off")
for ax,c in zip(axes,order):
    sub=H[H.cohort==c]; alluvial(ax,sub,["DRB345","HLA-DRB1","HLA-DQA1","HLA-DQB1"],"DRhap",f"{SHORT[c]} (n={len(sub)})")
fig.suptitle("Population-level class II haplotype flow (DRB3/4/5 - DRB1 - DQA1 - DQB1), coloured by secondary DRB gene",fontsize=12)
fig.tight_layout(); fig.savefig("figures/fig7_classII_flow_by_cohort.png",dpi=140)
feat=pd.get_dummies(H[["DRB345","HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DQA2","HLA-DQB2"]].astype(str)).astype(float)
X=feat.values-feat.values.mean(0); U,S,Vt=np.linalg.svd(X,full_matrices=False); pc=U[:,:2]*S[:2]; ev=S**2/np.sum(S**2)
H["pc1"],H["pc2"]=pc[:,0],pc[:,1]; H["sample"]=[h.split("#")[0] for h in H.index]
fig,(ax,ax2)=plt.subplots(1,2,figsize=(19,8.5)); cmap={"DRB3":"#1f77b4","DRB4":"#ff7f0e","DRB5":"#2ca02c","none":"#7f7f7f"}
for s,g in H.groupby("sample"):
    if len(g)==2: ax.plot(g.pc1,g.pc2,color="#bbbbbb",lw=0.5,zorder=1)
for k,col in cmap.items():
    g=H[H.DRhap==k]; ax.scatter(g.pc1,g.pc2,s=14,color=col,label=f"DRA-{k}-DRB1-DQA1" if k!="none" else "DRA-DRB1-DQA1 (no DRB3/4/5)",zorder=2,alpha=0.8)
for grp,g in H.groupby("DRB1grp"):
    if len(g)>=8: ax.text(g.pc1.median(),g.pc2.median(),grp,fontsize=8,ha="center",va="center",bbox=dict(boxstyle="round,pad=0.15",fc="white",ec="none",alpha=0.7))
ax.set_xlabel(f"PC1 ({ev[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({ev[1]*100:.1f}%)"); ax.legend(fontsize=8)
ax.set_title("Diplotype PCA on class II two-field alleles (DRB3/4/5, DRB1, DQA1, DQB1, DQA2, DQB2)\nlines join the two haplotypes of an individual; labels = DRB1 allele group",fontsize=10)
for c in present(H.cohort)[::-1]:
    g=H[H.cohort==c]; ax2.scatter(g.pc1,g.pc2,s=16 if c.startswith("HPRC-Rest") else 24,color=COLOR[c],label=f"{SHORT[c]} (n={len(g)})",alpha=0.55 if c=="HPRC-Rest" else 0.9,zorder=1 if c=="HPRC-Rest" else 2,edgecolors="none")
ax2.set_xlabel(ax.get_xlabel()); ax2.set_ylabel(ax.get_ylabel()); ax2.legend(fontsize=8); ax2.set_title("Same PCA coloured by cohort",fontsize=10)
fig.tight_layout(); fig.savefig("figures/fig8_classII_diplotype_pca.png",dpi=160)
u=H[(H.DRB1grp=="DRB1*13")&(H["HLA-DQA1"].str.startswith("DQA1*01"))&(H["HLA-DQB1"].str.startswith("DQB1*05"))]
print("DRB1*13-DQA1*01-DQB1*05 haplotypes:\n", u[["cohort","DRB345","HLA-DRB1","HLA-DQA1","HLA-DQB1"]].to_string())
