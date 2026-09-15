#!/usr/bin/env python3
"""Figure 4: candidate novel coding HLA alleles with independent-assembly and read support.
Run inside the hla-viz directory (data/novel_coding_support.tsv, data/rs/*.novel_pileup.tsv)."""
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sup=pd.read_csv("data/novel_coding_support.tsv",sep="\t")
POP=pd.read_csv("data/hprc_r2_populations.tsv",sep="\t").set_index("sample")
def hprc(s): return f"{POP.analysis_group[s]} ({POP.population[s]})"
# Immuannot cds_mut letters are assembly-then-known; the extraction script labelled them the other way round
sup=sup.rename(columns={"known":"assembly_base","asm":"known_base"})
sup["ok"]=sup.assembly_base==sup.asm_base_check
alleles=[("apr003","1","HLA-A","APR"),("apr011","1","HLA-C","APR"),("apr001","1","HLA-DRB1","APR"),("apr038","1","HLA-DRB1","APR"),("apr002","1","HLA-DQB1","APR"),("apr015","2","HLA-DPB1","APR"),("HG02717","1","HLA-DQB1",hprc("HG02717")),("NA20346","2","HLA-DPA1",hprc("NA20346")),("ksa004","1","HLA-A","JaSaPaGe-Saudi"),("HG03516","2","HLA-DPA1",hprc("HG03516"))]
indel_only={("ksa004","HLA-A"):("A*02:01:01:134Q","3 indels only (1-bp deletion, 2 insertions)"),("HG03516","HLA-DPA1"):("DPA1*03:01:01:01","1-bp insertion only")}
with_indels={("apr038","HLA-DRB1")}
reads={("HG02717","HLA-DQB1"):(8,14),("NA20346","HLA-DPA1"):(24,51)}
# checked against IPD-IMGT/HLA 3.65 (Immuannot uses 3.55): HG02717's coding sequence (exons 2-6) equals DQB1*02:180:02,
# created in release 3.56 (Feb 2024), and the HPRC truth set (Lai et al. 2024) calls DQB1*02:180; NA20346's DPA1 has no
# exact coding match in 3.65 (closest DPA1*03:02:02, 1 substitution in exon 2)
named_later={("HG02717","HLA-DQB1"):"Ala>Asp vs 02:02:01; = DQB1*02:180:02 (named in IPD-IMGT/HLA 3.56), not novel",
             ("NA20346","HLA-DPA1"):"Ala>Met vs 01:03:01; still novel in 3.65 (closest 03:02:02, 1 subst.)"}
rows=[]
for s,h,g,c in alleles:
    x=sup[(sup["sample"]==s)&(sup.gene==g)]
    if (s,g) in indel_only:
        cl,desc=indel_only[(s,g)]; rows.append((s,h,g,c,cl,desc,"-","-","no public reads" if not c.startswith("HPRC") else "-",True)); continue
    x=x[x.ok]
    priv=(x.n_other_haps_with_31mer==0).sum(); oth=x.other_hap_of_same_individual.all()
    import re as _re
    aa=", ".join(_re.sub(r"\((\w+)\)","",a.split(",")[0]).replace("Rrg","Arg").replace("Tre","Thr").replace("<",">") for a in x.aa)
    desc=f"{len(x)} subst."+(" + indels" if (s,g) in with_indels else "")+f": {aa}"
    rd=reads.get((s,g)); rdtxt=f"{rd[0]}/{rd[1]} assembly base, {rd[1]-rd[0]} other hap" if rd else "no public reads"
    if (s,g) in named_later: desc=named_later[(s,g)]
    rows.append((s,h,g,c,x.closest.iloc[0].replace("HLA-",""),desc,str(priv),"yes" if oth else "no",rdtxt,(s,g) in with_indels))
tab=pd.DataFrame(rows,columns=["sample","hap","gene","cohort","closest","desc","private","other hap","reads","indel"])
tab.drop(columns="indel").to_csv("data/novel_alleles_table.tsv",sep="\t",index=False)
fig=plt.figure(figsize=(17,11)); gs=fig.add_gridspec(2,2,height_ratios=[0.9,1],hspace=0.45,wspace=0.15)
ax=fig.add_subplot(gs[0,:]); ax.axis("off")
cell=[[r["sample"]+"#"+str(r["hap"]),r["cohort"],r["gene"].replace("HLA-",""),r["closest"],r["desc"][:100],r["private"],r["other hap"],r["reads"]] for _,r in tab.iterrows()]
t=ax.table(cellText=cell,colLabels=["haplotype","cohort","gene","closest known allele","coding difference(s) (known > assembly)","private\nvariants","also on\nother hap","own 1000G Illumina reads\n(realigned to assembly)"],loc="upper center",cellLoc="left",colWidths=[0.08,0.1,0.05,0.13,0.34,0.06,0.06,0.18],bbox=[0,0.05,1,0.9])
t.auto_set_font_size(False); t.set_fontsize(7.5)
for (i,j),c in t.get_celld().items():
    if i==0: c.set_text_props(weight="bold"); c.set_facecolor("#e8e8e8")
    elif tab.indel.iloc[i-1]: c.set_facecolor("#fde0dd")
    elif tab.reads.iloc[i-1][0].isdigit(): c.set_facecolor("#e0f3db")
ax.set_title("a  Candidate novel coding alleles (Immuannot ':new' with CDS differences) in classical HLA genes.  private = variant 31-mer absent from all other 609 haplotypes;\n    red = differences include indels (typical HiFi homopolymer assembly errors, treat as artefacts); green = heterozygous support in the individual's own short reads",loc="left",fontsize=9.5)
def pileup(ax,fn,pos,title):
    p=pd.read_csv(fn,sep="\t",header=None,names=["pos","ref","depth","refn","A","C","G","T"])
    p=p[(p.pos>=pos-40)&(p.pos<=pos+40)].reset_index(drop=True)
    x=np.arange(len(p)); bottom=np.zeros(len(p)); cols={"A":"#2ca02c","C":"#1f77b4","G":"#ff7f0e","T":"#d62728"}
    for b in "ACGT":
        v=np.where(p.ref==b,p.refn,0)+p[b].values
        ax.bar(x,v,bottom=bottom,color=cols[b],width=0.9,lw=0); bottom+=v
    i=int(np.where(p.pos.values==pos)[0][0]); ax.annotate("novel codon position",xy=(i,bottom[i]),xytext=(i,bottom.max()*1.15),ha="center",fontsize=8,arrowprops=dict(arrowstyle="->",color="k"))
    ax.set_xticks(x[::10]); ax.set_xticklabels([str(v) for v in p.pos.values[::10]],fontsize=7,rotation=45); ax.set_ylabel("reads"); ax.set_title(title,loc="left",fontsize=9.5); ax.set_ylim(0,bottom.max()*1.3)
    ax.set_xlabel("position on the assembled MHC contig (bp); mixed columns = heterozygous sites",fontsize=8)
    for j in range(len(p)):
        c=p.iloc[j]; alt={b:c[b] for b in "ACGT" if b!=c.ref}
        a,n=max(alt.items(),key=lambda kv:kv[1])
        if n>=3 and n/c.depth>=0.2: ax.text(j,bottom[j]+1,f"{c.ref}/{a}",fontsize=6,ha="center",rotation=90)
    ax.legend(handles=[plt.Rectangle((0,0),1,1,color=cols[b]) for b in "ACGT"],labels=list("ACGT"),fontsize=7,loc="upper right",ncol=4)
pileup(fig.add_subplot(gs[1,0]),"data/rs/HG02717.novel_pileup.tsv",4236512,"b  HG02717 hap1 HLA-DQB1, Ala>Asp at CDS 266 (gene on - strand): the allele is DQB1*02:180:02 (IPD-IMGT/HLA 3.56)\n    8/14 reads carry the assembly base, 6/14 the other haplotype")
pileup(fig.add_subplot(gs[1,1]),"data/rs/NA20346.novel_pileup.tsv",4622149,"c  NA20346 hap2 HLA-DPA1, Ala>Met at CDS 124 (gene on - strand): novel in IPD-IMGT/HLA 3.65\n    24/51 reads carry the assembly base, 27/51 the other haplotype")
fig.savefig("figures/fig4_novel_coding_alleles.png",dpi=160,bbox_inches="tight")
print(tab.drop(columns="indel").to_string())
