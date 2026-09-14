# Literature background

## Ito-Naito et al., HLA 2026;108(2):e70898 (PMC13457674)

*Genetic Diversity and Haplotype Architecture of Non-Classical HLA Class I Genes
in the Japanese Population.* Full-length HLA-E (~8.4 kb), HLA-F (~4.8 kb) and
HLA-G (~4.9 kb) amplicons (promoter to 3' UTR) sequenced with PacBio in 531
Japanese individuals (1,062 alleles); four-field typing with LAA/pbAA +
NGSengine, TG-repeat genotyping with TRGT, LD with PLINK, frequencies with
PyPop/BIGDAWG.

Findings that matter for this workflow:

- 36 novel and 16 extended alleles, nearly all in noncoding/regulatory
  sequence; protein-level diversity is tiny (F*01:01 96.7 %, E*01:03 69.2 %,
  G*01:01 53.1 % / G*01:04 45.4 %). Assemblies cover full genes, so our calls
  should be reported at four-field resolution and novel alleles split into
  coding vs noncoding (`allele_4field`, `novel_class` in `hla_calls.tsv`,
  `plots/novel_allele_class.png`).
- HLA-G 14-bp 3' UTR indel (rs371194629) co-segregates with four-field
  alleles (deletion with G*01:04); rs2523405-T with F*01:01:02 / F*01:01:06.
  Worth checking in the per-gene pgr-tk bundle plots (3' end of HLA-G).
- No LD between HLA-E and F/G, moderate F–G LD (F*01:01:02:09 ~ G*01:01:03:03,
  D' = 1.0). Haplotype-resolved assemblies give F–G phase directly.
- Homopolymer tracts in HLA-E intron 4 cause "repeat region ambiguities";
  assembly-based calls need the same caution (Immuannot `template_distance`
  small but nonzero).
- Reference sets for East Asian and other under-represented populations are
  the stated gap; the APR, JaSaPaGe and HPRC r2 cohorts here address it.
