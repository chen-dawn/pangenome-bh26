#!/usr/bin/env python3
"""Figure 3: MHC homozygosity / collapsed haplotypes. Inputs (hla-viz/data): het_density_100kb.tsv (1000G phased SNPs),
rs/<sample>.<hap>.sites.tsv (own 1000G reads realigned to own assembled haplotype 1, scripts/read_support.sh),
rs/<sample>.h1h2.paf (minimap2 asm20 hap1 vs hap2 with cs tags)."""
import pandas as pd, numpy as np, re, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
dens=pd.read_csv("data/het_density_100kb.tsv",sep="\t",index_col=0)
samples=[("NA18976","HPRC r2 NA18976 (JPT)"),("NA19909","HPRC r2 NA19909 (ASW)"),("HG00544","HPRC r2 HG00544 (CHS)"),("NA18952","NA18952 (JPT)"),("NA18940","NA18940 (JPT)")]
genes38={"HLA-A":29942532,"HLA-C":31268749,"HLA-B":31353872,"C4A":31982057,"HLA-DRB1":32578775,"HLA-DQB1":32659467,"HLA-DPB1":33075990}
def hetsites(fn,lo=0.25,hi=0.75,mind=15):
    s=pd.read_csv(fn,sep="\t",header=None,names=["ctg","pos","ref","depth","refn","alt","altn","frac"])
    return s[(s.depth>=mind)&(s.frac>=lo)&(s.frac<=hi)]
tracks=[("NA18976","data/rs/NA18976.NA18976_1.sites.tsv","HPRC NA18976 hap1"),("NA19909","data/rs/NA19909.NA19909_1.sites.tsv","HPRC NA19909 hap1"),("HG00544","data/rs/HG00544.HG00544_1.sites.tsv","HPRC HG00544 hap1"),("NA18952","data/rs/NA18952.NA18952_1.sites.tsv","HPRC NA18952 hap1"),("NA18952j","data/rs/NA18952j.NA18952.JaSaPaGe_1.sites.tsv","JaSaPaGe NA18952 hap1"),("NA18940","data/rs/NA18940.NA18940_1.sites.tsv","HPRC NA18940 hap1")]
def paf_div(fn):
    win={}
    for line in open(fn):
        f=line.rstrip().split("\t")
        if int(f[11])<30 or int(f[10])<20000: continue
        cs=[x for x in f if x.startswith("cs:Z:")]
        if not cs: continue
        tpos=int(f[7])
        for m in re.finditer(r"(:\d+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+)",cs[0][5:]):
            t=m.group(0)
            if t[0]==":": tpos+=int(t[1:])
            elif t[0]=="*": win[tpos//50000]=win.get(tpos//50000,0)+1; tpos+=1
            elif t[0]=="-": tpos+=len(t)-1
    return win
cols={"NA18976":"#d62728","NA19909":"#ff7f0e","HG00544":"#9467bd","NA18952":"#1f77b4","NA18952j":"#17becf","NA18952.JaSaPaGe":"#17becf","NA18940":"#7f7f7f"}
fig,axes=plt.subplots(3,1,figsize=(15,11),gridspec_kw={"height_ratios":[1.2,1.2,1]})
ax=axes[0]
for s,lab in samples: ax.plot(dens.index/1e6,dens[s]+0.5,lw=1.3,label=lab,color=cols[s])
ax.set_yscale("log"); ax.set_ylabel("heterozygous SNPs per 100 kb\n(1000G Illumina, GRCh38)"); ax.set_xlim(28.5,33.5)
for g,p in genes38.items(): ax.axvline(p/1e6,color="k",lw=0.4,ls=":"); ax.text(p/1e6,3000,g.replace("HLA-",""),fontsize=7,ha="center")
ax.legend(fontsize=8,ncol=3,loc="lower right"); ax.set_title("a  Short-read heterozygosity across the MHC (1000G phased genotypes, independent of the assemblies)",loc="left",fontsize=11); ax.set_xlabel("GRCh38 chr6 (Mb)")
ax=axes[1]
for s,fn,lab in tracks:
    h=hetsites(fn); w=(h.pos//50000)*50000; c=h.groupby(w).size()
    ax.plot(c.index/1e6,c.values+0.5,lw=1.3,label=lab,color=cols[s],ls="--" if "JaSaPaGe" in lab else "-")
ax.set_yscale("log"); ax.set_ylabel("sites with 25-75% alternative base\nper 50 kb"); ax.set_xlabel("assembled MHC haplotype 1 (Mb from start of extracted region)")
ax.legend(fontsize=8,ncol=3,loc="lower right"); ax.set_title("b  1000G reads of the same individual realigned to its own assembled haplotype 1 (pileup heterozygosity)",loc="left",fontsize=11)
ax=axes[2]
for s,lab in [("NA18976","HPRC NA18976"),("NA19909","HPRC NA19909"),("HG00544","HPRC HG00544"),("NA18952","HPRC NA18952"),("NA18952.JaSaPaGe","JaSaPaGe NA18952"),("NA18940","HPRC NA18940")]:
    w=paf_div(f"data/rs/{s}.h1h2.paf"); k=sorted(w)
    ax.plot([x*50000/1e6 for x in k],[w[x]+0.5 for x in k],lw=1.3,label=lab,color=cols[s],ls="--" if "JaSaPaGe" in lab else "-")
ax.set_yscale("log"); ax.set_ylabel("substitutions hap1 vs hap2\nper 50 kb (minimap2 asm20)"); ax.set_xlabel("assembled MHC haplotype 1 (Mb)"); ax.legend(fontsize=8,ncol=3,loc="lower right")
ax.set_title("c  Divergence between the two assembled haplotypes of each individual",loc="left",fontsize=11)
fig.tight_layout(); fig.savefig("figures/fig3_mhc_homozygosity.png",dpi=160)
