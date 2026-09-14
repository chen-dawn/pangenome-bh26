# K-PanRef: Korean pangenome reference (14 individuals)

Paper: Shin et al., "A Korean pangenome reference of 14 healthy individuals supports structural variant analysis in disease genomes",
       medRxiv 2026, doi:10.64898/2026.07.06.26357367
Code: https://github.com/Da-Blessed/K-PanRef
Graph: Minigraph-Cactus v3.0.0, references CHM13 + GRCh38.
Downloaded: 2026-09-14 (biohackathon) by Robert Hoehndorf

- KPanRef.gbz     graph (GBZ, vg giraffe ready)   } Zenodo https://doi.org/10.5281/zenodo.20810335
- KPanRef.gfa.gz  graph (GFA)                     }
- KPanRef.vcf.gz  graph VCF                       }

Assemblies: deposited in K-BDS KNA under KAP242400 (reads: KRA KAP242397), but as of 2026-09-14
the accession is not searchable/released in KNA, so assembly FASTAs are NOT here.
Haplotype paths can be extracted from the GBZ (vg paths -x KPanRef.gbz -S <sample> -F) as a fallback.
