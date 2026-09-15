# HLA pangenome figures (754 haplotypes)

Cohorts: APR (UAE Arab) 106; JaSaPaGe Saudi 18 and Japanese (1000G JPT) 20; HPRC release 2 split by population into
Japanese (JPT) 32, Ashkenazi Jewish 2 (HG002 only), other East Asian 70 and rest 360; K-PanRef Korean 28 and CPC Chinese 116
(haplotype sequences read from their pangenome graphs); GRCh38 and CHM13. All calls are Immuannot (IPD-IMGT/HLA 3.55) on
MHC regions cut with pgr-tk; per-gene sequences are cut at the Immuannot coordinates; bundle analyses use pgr-tk 0.5.1.
Scripts: `hla/analysis/` in `leechuck/pangenome-bh26` (run on the NIG cluster with `hla/nig/make_figures.sh`); cohort
order and colors in `analysis/cohorts.py`; data tables in `tables/`.

**fig2_mhc_extraction.png - The MHC comes out complete from 747 of 752 haplotypes.**
Left: fraction of the GRCh38 MHC (chr6:28,510,120-33,480,577) covered by the extracted segments of each haplotype; numbers
give haplotypes with coverage >= 0.99 per cohort. Right: number of MHC segments per haplotype (720 of 752 on one segment).
Below 0.99: four CPC Chinese haplotypes (lowest 0.92) and one JaSaPaGe Saudi haplotype (0.96).

**fig3_mhc_homozygosity.png - Two HPRC assemblies are genuinely MHC-homozygous; one JaSaPaGe assembly is a duplicated haplotype.**
(a) Heterozygous SNPs per 100 kb across the MHC from 1000G high-coverage Illumina genotypes (assembly-independent).
NA18976 (JPT) and NA19909 (ASW) carry <100 heterozygous sites outside the paralog-rich HLA-A/-H and DRB windows,
versus >10,000 for typical individuals. (b) The same individuals' 1000G reads realigned to their own assembled
haplotype 1: sites with 25-75% alternative base per 50 kb. Both haplotypes of JaSaPaGe NA18952 (dashed) show the
full heterozygosity of the individual, i.e. the second haplotype is missing from that assembly, whereas the HPRC
assembly of the same individual separates the two. (c) Substitutions between the two assembled haplotypes of each
individual (minimap2 asm20): NA18976 202 and NA19909 38 differences over 5.2 Mb; JaSaPaGe NA18952 hap1 = hap2 (10 differences).

**fig4_novel_coding_alleles.png - Candidate novel coding alleles and their support.**
(a) Immuannot ":new" alleles (absent from IPD-IMGT/HLA 3.55) with CDS differences in classical genes, from the
610-haplotype run. "private" = number of substitutions whose 31-mer occurs in no other of those 610 haplotypes; most
substitutions of the multi-change alleles (apr003 A*24, apr011 C*03, apr001 DRB1*16) are seen in hundreds of other
haplotypes, i.e. these look like recombinant alleles of known segments rather than errors. Rows in red differ from the
closest allele only by 1-bp indels (HiFi homopolymer artefacts). (b, c) Base-level pileups of the individual's own 1000G
Illumina reads realigned to the assembled contig around the codon: HG02717 DQB1 Ala>Asp (8/14 reads carry the assembly
base, the rest the other haplotype) and NA20346 DPA1 Ala>Met (24/51). Checked against IPD-IMGT/HLA 3.65: the HG02717 coding
sequence (exons 2-6) equals DQB1*02:180:02, created in release 3.56 (February 2024; the HPRC truth set of Lai et al. 2024
also calls DQB1*02:180), so it is not novel; NA20346 DPA1 remains novel (closest DPA1*03:02:02, one substitution in exon 2).
APR/JaSaPaGe reads are not public, so those rows are assembly-only.

**fig5_population_hla.png - HLA allele landscape by cohort.**
Top: two-field allele frequencies (top 6 per cohort) for HLA-A, -B, -DRB1; e.g. A*24:02 in 35% (JaSaPaGe) and 41% (HPRC)
of Japanese haplotypes, A*11:01 in 29% of CPC Chinese and 26% of other East Asian, B*51:01 in 28% of Saudi, DRB1*03:01 in
18% of APR, DRB1*12:02 in 16% of CPC Chinese. Bottom: secondary DRB gene (DR haplogroup) frequencies, C4A/C4B long/short
forms, and the fraction of gene copies whose full-length sequence is absent from IPD-IMGT/HLA 3.55 (57-79% for DRB1).
The Jewish column is a single individual (HG002, 2 haplotypes) and is not a population estimate.

**fig6_classII_haplotype_flow.png / fig7_classII_flow_by_cohort.png - Gene-level class II haplotype flow (after Chin, ASHI 2023).**
Alluvial plot of two-field alleles along DRB3/4/5 - DRB1 - DQA1 - DQB1 - DQA2 - DQB2 - TAP2 - TAP1 for all 752
haplotypes (ribbons colored by secondary DRB gene), and per cohort for the DR-DQ block. 518 distinct gene-level
strings among 752 haplotypes; the DRB1*13-DQA1*01-DQB1*05 combination Chin flagged in HG03516/NA18906 recurs in 9
HPRC haplotypes. DRB4 haplotypes make up 44-45% of Saudi and JaSaPaGe Japanese and 39% of Korean haplotypes; DRB3 46% of
CPC Chinese.

**fig8_classII_diplotype_pca.png - Diplotype PCA on class II two-field alleles.** Each point is a haplotype, grey lines
join the two haplotypes of an individual; left: color = secondary DRB gene, labels = DRB1 allele group; right: cohort.

**fig9_classII_dendrogram.png - pgr-tk bundle-distance dendrogram of the DRA..DMA region (753 haplotypes).**
Color strips: DRB1 allele group, secondary DRB gene, cohort. Clusters follow DRB1 allele groups; within-cluster share of
the most common 2-field DRB1 allele is 39% at 10 clusters and 52% at 80 (`tables/classII_cluster_purity.tsv`).

**fig10_classII_bundle_graph.png - Principal-bundle graph of the class II region.** 890 bundles from
pgr-pbundle-decomp (w=48, k=56, r=2, min_span=8); nodes = bundles (size = length, color = number of haplotypes),
edges = consecutive bundles on a haplotype. The DRB block forms the tangled part; DQ/DO/TAP/DM is nearly linear.

**fig11_classII_bundle_pca.png - Diplotype PCA on principal-bundle presence (sequence level).** The 751 x 890 bundle
presence matrix separates DRB4 (DR7/DR9/DR4), DRB5 (DR15/DR16), DRB3 (DR3/11/12/13/14) and DRB1-only (DR1/DR8/DR10)
haplogroups on PC1/PC2 (27% and 16% of variance); right panel: the same PCA colored by cohort, showing that the
haplogroup clusters contain all cohorts. Grey lines join the two haplotypes of each individual.

**fig12_typing_concordance.png - Typing methods agree at 2-field resolution.**
Fraction of individuals with identical unordered genotypes, per gene, at 1-4 fields. Left: FuFiHLA (IPD-IMGT/HLA 3.65, run
on tiled pseudo-reads of both assembled haplotypes) vs Immuannot on the same assemblies, 373-376 individuals; 2-field
agreement 92% (DRB1) to 100%. Of the DRB1 disagreements, 19 are alleles both methods call novel and 10 are FuFiHLA
artefacts (a second allele invented on a homozygous gene, absent from the assembly). Middle and right: T1K on 1000G
Illumina MHC reads vs Immuannot and vs FuFiHLA for the 23 assembled individuals typed by T1K so far. Not plotted: T1K
agrees with Sanger typing of 412 1000G samples (Gourraud et al. 2014) in 98.8-100% of genes at 2 fields, and against the
HPRC 4-field labels of Lai et al. 2024 (44 individuals) Immuannot reaches >= 97.7% and FuFiHLA 100% at 2 fields
(`tables/typing_concordance*.tsv`).

**Interactive / per-gene:** `results/<GENE>.html` (pgr-tk bundle plots), `classII/classII.html` (DRA..DMA),
`results/MHC.html` (whole MHC), `results/<GENE>.viz.png` (pggb + odgi), `results/*.png` (Immuannot summaries).
