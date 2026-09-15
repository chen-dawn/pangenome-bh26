# HLA gold-standard datasets

This folder holds **experimental / high-confidence HLA gold-standard labels
only** — datasets backed by wet-lab typing or targeted-capture-verified
assemblies, not computational/in-silico HLA calls. The earlier Abi-Rached et
al. 2018 dataset (2,693 1000G samples, exome-based PolyPheMe *computational*
typing) is deliberately excluded for this reason.

## Quick answer

| Dataset | Folder | Samples | Genes | Resolution | What the labels represent |
|---|---|---:|---|---|---|
| **Gourraud et al. 2014 / 1000 Genomes** | `1000G_2014/` | 1,267 | A, B, C, DRB1, DQB1 | ARS-exon based (exons 2/3 class I, exon 2 class II); ambiguous calls preserved | **Experimental Sanger HLA typing** |
| **Lai et al. 2024 / HPRC** | `HPRC_2024/` | 44 | A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1, DRB3, DRB4, DRB5 | **4-field** | High-confidence labels from phased HPRC assemblies, cross-checked with targeted HLA capture sequencing |

Use **Gourraud 2014** for a conventional experimental benchmark. Use
**Lai 2024 / HPRC** when 4-field labels are required (e.g. to check Immuannot
`consensus` calls at full resolution instead of truncating to 2-field).

---

## 1. `1000G_2014/` — Gourraud et al. 2014, experimental Sanger HLA typing

- **Paper:** Gourraud PA, Khankhanian P, Cereb N, Yang SY, Feolo M, Maiers M,
  Rioux JD, Hauser S, Oksenberg J. "HLA Diversity in the 1000 Genomes
  Dataset." *PLOS ONE* 2014.
  DOI: [10.1371/journal.pone.0097282](https://doi.org/10.1371/journal.pone.0097282) ·
  [PMC4079705](https://pmc.ncbi.nlm.nih.gov/articles/PMC4079705/) ·
  [IGSR FAQ](https://www.internationalgenome.org/faq/was-hla-diversity-studied-in-igsr/)
- **Samples:** 1,267 individuals, 14 populations, 4 major ancestral groups.
  Sample IDs are standard 1000G/Coriell IDs (`NA*`).
- **Genes:** HLA-A, -B, -C, -DRB1, -DQB1.
- **Method:** Locus-/group-specific PCR + **Sanger sequencing** (SBT),
  covering exons 2 and 3 for class I and exon 2 for class II only — real
  wet-lab typing, not full-gene sequence truth. Compared against IMGT
  2.26.0. Typing ambiguities (identical ARS-exon sequences across alleles)
  are preserved in the data rather than force-resolved, e.g.
  `32:01:01/32:01:02`.
- **Downloaded from:**
  `https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/`
- **Files on disk:**
  - `20140702_hla_diversity.txt` (371,272 bytes) — genotype table, columns
    `id, sbgroup, A, A.1, B, B.1, C, C.1, DRB1, DRB1.1, DQB1, DQB1.1`
    (`id` = Coriell/1000G sample ID, `sbgroup` = sub-population, paired
    columns = the two alleles per locus, ambiguous calls kept as `/`-joined
    strings).
  - `README_20140702_hla_diversity.txt` (3,870 bytes) — original README with
    full methods.
- **Limitation:** resolution is capped at what the sequenced ARS exons
  support — do not upgrade an ambiguous 2014 call to a modern 3-/4-field
  allele by inference.

## 2. `HPRC_2024/` — Lai et al. 2024, 4-field HPRC gold-standard labels

- **Paper:** Lai SK, Luo AC, Chiu IH, Chuang HW, Chou TH, Hung TK, Hsu JS,
  Chen CY, Yang WS, Yang YC, Chen PL. "A novel framework for human leukocyte
  antigen (HLA) genotyping using probe capture-based targeted next-generation
  sequencing and computational analysis." *Comput Struct Biotechnol J* 2024.
  DOI: [10.1016/j.csbj.2024.03.030](https://doi.org/10.1016/j.csbj.2024.03.030) ·
  [PMC11035020](https://pmc.ncbi.nlm.nih.gov/articles/PMC11035020/)
- **Samples:** 44 HPRC samples (IDs like `HG002`, `HG00438`, ...).
- **Genes:** HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1, -DRB3, -DRB4,
  -DRB5 — a superset of Immuannot's class I/II gene set.
- **Method:** HLA annotated on phased HPRC personal assemblies by two
  independent approaches, discrepancies manually inspected on the assembly
  contigs, and cross-checked against targeted HLA capture sequencing on DNA
  from the same individuals. Supplementary-6 is the resulting label table
  (not the OptiType/SpecHLA/HLA-VBSeq prediction outputs elsewhere in the
  supplement — those are methods being *evaluated* in the paper, not truth).
- **Downloaded from:** Europe PMC's supplementary-files bundle for
  PMC11035020 (`https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11035020/supplementaryFiles`),
  which mirrors the publisher's Supplementary-6 without the NCBI PMC
  download gate. Verified byte-for-byte size match (16,951 bytes) against
  the file NCBI PMC lists for `mmc6.xlsx`.
- **File on disk:** `Lai2024_HPRC_Supplementary-6_HLA_4field_genotypes.xlsx`
  — one row per sample (`Sample` column, e.g. `HG002`), one column per gene
  (`HLA-A`, `HLA-B`, ..., `HLA-DRB5`), 4-field alleles per cell (e.g.
  `A*26:01:01:01`).

---

## Using this against the pangenome/Immuannot pipeline

- Both sets use standard 1000G/HPRC sample IDs (`NA*`/`HG*`), so they line up
  directly with samples in `hla/nig/inputs-nig.yml` cohorts that trace back
  to 1000G/HPRC individuals (e.g. `HPRC_r2`, and the JPT trios shared between
  `HPRC_r2` and `JaSaPaGe`).
- `HPRC_2024/` gives 4-field truth for the *same assemblies this pipeline
  already processes* (HPRC_r2) — it's the closer match for validating
  Immuannot's `consensus` calls in `hla_calls.tsv` at full resolution without
  truncation.
- `1000G_2014/` is exon-only Sanger typing, independent of any assembly
  method, but capped below 4-field resolution — truncate Immuannot calls to
  match (e.g. `A*02:01:01:01` → `02:01`) and expect some entries to only be
  comparable at the ambiguous multi-allele level reported in the file.
