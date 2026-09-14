# Arab Pangenome Reference (APR / UPR)

Source: https://www.mbru.ac.ae/the-arab-pangenome-reference/
Paper: Nassir et al., "A draft UAE-based Arab pangenome reference", Nat Commun 2025, doi:10.1038/s41467-025-61645-w
Code: https://github.com/muddinmbru/arab_pangenome_reference
Downloaded: 2026-09-14 (biohackathon) by Robert Hoehndorf

- assemblies/  106 haplotype-resolved polished assemblies (53 individuals x 2 haplotypes), *.polished.fa.gz.
               From SharePoint folder "APR Assemblies" (anonymous share link on the MBRU page).
- graph/       Minigraph-Cactus pangenome graph (CHM13-based, "apr_review_v1_2902_chm13"):
               .gbz .gfa.gz .hapl .d9.dist .d9.snarls .min .vcf.gz
               From SharePoint folder "APR Nuclear/Pangenome".
- set02.zip    Zenodo record https://zenodo.org/records/13752609 (42.8 GB): contains 53 of the same
               *.polished.fa.gz assemblies (a subset of assemblies/, i.e. redundant). See set02.zip.listing.txt.

Not downloaded: raw reads, per-sample GRCh38 DeepVariant VCFs ("APR Nuclear/Variants"), mitochondrial folder, annotation folder.
