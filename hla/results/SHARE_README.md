# HLA pangenome analysis, BioHackathon Japan 2026 (shared folder `/home/asianhla/data/upload/HLA/`)

Generated 14-15 September 2026 by Robert Hoehndorf (KAUST) with Claude Code on the NIG supercomputer
(partition `asianhla-c32`, node asianhla-vm, 28 cores / 200 GB). Everything here was produced by the CWL
workflow and scripts in `workflow/` (git: https://github.com/leechuck/pangenome-bh26, directory `hla/`).
Questions: robert.hoehndorf@kaust.edu.sa

## Inputs

610 haplotype assemblies (`inputs-610-haplotypes.yml` lists every file with sample, haplotype and cohort):

| cohort | haplotypes | source on this cluster |
|---|---|---|
| APR (UAE Arab Pangenome Reference) | 106 | `/home/asianhla/data/upload/APR/assemblies` |
| HPRC_r2 | 464 | `/home/asianhla/data/HPRC_r2/fasta` |
| JaSaPaGe-Saudi (ksa001-009) | 18 | `/home/asianhla/data/JaSaPaGe/assembly_clean/fasta` |
| JaSaPaGe-Japanese (1000G JPT NA*) | 20 | same folder |
| KPanRef-Korean (K-PanRef, 14 individuals) | 28 | graph only: `/home/asianhla/data/upload/KPanRef/KPanRef.gbz` |
| CPC-Chinese (CPC Phase 1, 58 individuals) | 116 | graph only: `/home/asianhla/data/CPC/CPC.Phase1.CHM13v2-full.gfa` |
| REF (GRCh38, CHM13) | 2 | `/home/asianhla/data/HPRC_r2/fasta` |

HPRC r2 is analysed per population (`workflow/data/hprc_r2_populations.tsv`, from the 1000G sequence indexes and the HPRC
release 2 sample metadata): HPRC-Japanese (JPT, 32 haplotypes), HPRC-Jewish (HG002, Ashkenazi, the only Jewish individual
in HPRC r2: 2 haplotypes, not a population estimate), HPRC-EastAsian (CHB/CHS/CDX/KHV and HG005, 70), HPRC-Rest (360).

K-PanRef and CPC assemblies are not available as FASTA, so their MHC sequences come from the pangenome graphs:
`workflow/nig/kpanref_mhc.sbatch` (all haplotype paths of each sample from the GBZ with `vg paths -F`, segments with an
alignment of >= 50 kb matching bases, MAPQ >= 20, to the GRCh38 MHC) and `workflow/nig/cpc_mhc.sbatch` (odgi extract of
the CHM13 chr6 28-34 Mb reference nodes, then the span of every CPC haplotype walk through them read from the GFA, same
filter). HG00438/HG00621/HG00673 in the CPC graph are HPRC samples and are not added again. Both then go through the same
extraction, Immuannot and downstream steps (`workflow/nig/add_graph_cohorts.sbatch`).

Five JPT individuals (NA18940, NA18943, NA18945, NA18952, NA18970) are in both HPRC r2 and JaSaPaGe; the JaSaPaGe
copies are labelled `<sample>.JaSaPaGe`. **JaSaPaGe NA18952 is defective: both haplotype files carry the same MHC
haplotype** (the hap2 file even carries the HPRC contig name JBHIJT010000011.1); use the HPRC assembly instead.

## HLA typing beyond Immuannot

- `fufihla/` - FuFiHLA (Bioinformatics 2026, btag231) full-field typing of HLA-A/B/C/DRB1/DQA1/DQB1 for every assembled
  individual, from tiled pseudo-HiFi reads of both assembled haplotypes (README inside)
- `/home/asianhla/data/upload/1000G_MHC/` - all reads aligned to the MHC (same GRCh38 coordinates), chr6 alt haplotypes
  and HLA allele contigs for the 2504 1000 Genomes high-coverage samples, with T1K typing and the published 1000G HLA
  genotypes (README inside)
- `tables/typing_concordance*.tsv`, `figures/fig12_typing_concordance.png` - Immuannot vs FuFiHLA vs T1K vs published

## What is where

- `results/` - complete output of the workflow run "stage 4" (`hla_pangenome.cwl`, Toil, job 20565053), per haplotype
  and per gene. Key files:
  - `<sample>_<hap>.mhc.fa` / `.mhc.tsv`: extracted MHC segment(s) (PanSN names `sample#hap#contig:start-end[_rc]`) and extraction stats;
    `MHC.fa` = all 610 concatenated (3 Gb)
  - `<sample>_<hap>.gtf.gz`: Immuannot annotation (IPD-IMGT/HLA 3.55, IPD-KIR 2.13, RefSeq C4)
  - `hla_calls.tsv` (one row per annotated gene copy; use column `consensus`, `:new` = full-length sequence absent from
    IPD-IMGT/HLA; `novel_class` = known / novel_coding / novel_noncoding), `hla_calls_matrix.tsv`, `gene_copy_number.tsv`,
    `mhc_extraction_summary.tsv`, `allele_freq_<GENE>.png` and other summary plots
  - per gene (19 genes: HLA-A/B/C/E/F/G/DRA/DRB1/DRB3/DRB4/DRB5/DQA1/DQB1/DPA1/DPB1, MICA, MICB, TAP1, TAP2):
    `<GENE>.fa` (every annotated copy, gene orientation, 2 kb flank), `<GENE>.regions.tsv`, pgr-tk bundle outputs
    `<GENE>.bed/.ctg.summary.tsv/.mapg.gfa/.pmapg.gfa/.dist/.nwk/.ord/.svg/.html`, pggb graph `<GENE>.gfa/.og/.stats.tsv`,
    `<GENE>.viz.png`, `<GENE>.viz_depth.png`, `<GENE>.draw.png`, `pggb_<GENE>/` (pggb working directory)
  - whole MHC: `MHC.bed/.ctg.summary.tsv/.mapg.gfa/.pmapg.gfa/.dist/.nwk/.ord/.svg/.html` (pgr-tk bundles, w=80 k=56 r=6)
- `mhc_fastas/`, `immuannot_gtfs/` - the same MHC fastas and GTFs as flat folders (inputs for re-runs that skip
  extraction and annotation: `workflow/nig/stage3.sh <mhc_dir> <gtf_dir> <outname>`)
- `classII/` - DRA..DMA span from every haplotype (`classII.fa`, 461-650 kb each) and its pgr-tk bundle decomposition
  (`pgr-pbundle-decomp -w 48 -k 56 -r 2 --min-span 8 --bundle-length-cutoff 500 --bundle-merge-distance 2000`):
  `classII.svg/.html` (bundle plot), `.nwk` (dendrogram), `.ord` (bundle presence vectors), `.pmapg.gfa`, `.bed`
- `read_support/` - 1000G high-coverage Illumina reads of NA18976, NA19909, HG00544, NA18940, NA18952, HG02717, NA20346
  (MHC region only, `<sample>.mhc.bam`, from the public EBI CRAMs, `fetch.sh`/`fetch2.sh`), realigned to the individual's
  own assembled haplotype(s) (`run.sh`; `<sample>.<hap>.bam`), per-site allele counts (`*.sites.tsv`), depth per 50 kb,
  hap1-vs-hap2 alignments (`*.h1h2.paf`), base-level pileups at novel codons (`*.novel_pileup.tsv`), 1000G phased MHC
  genotypes of 9 samples (`mhc_9samples.vcf.gz`), Immuannot coding differences of the novel alleles (`novel_coding.tsv`)
- `figures/` - final figures with `CAPTIONS.md`; `tables/` - novel-allele support, class II haplotype strings, cluster purity,
  1000G heterozygosity tracks; `slides/` - summary deck (`hla_summary.pdf`, results and methods)
- `workflow/` - copy of the repository directory `hla/` (CWL workflow, tools, scripts, analysis scripts, README)

## How it was generated

1. Environment (`workflow/nig/setup_env.sh`): micromamba env with minimap2 2.28, samtools 1.21, pggb 0.7.4, odgi 0.9.2,
   seqkit, bcftools; pgr-tk 0.5.1 release binaries; Immuannot v3 + Data-2024Feb02; Toil 9.5 (cwltool's `--parallel`
   deadlocks on 610-way scatters). No containers on the node, so `--no-container`.
2. `workflow/scripts/make_inputs.py` builds the job file from the assembly folders (naming rules in the script).
3. `sbatch workflow/nig/run_nig.sbatch <job.yml> <outdir>` runs `workflow/hla_pangenome.cwl` (steps and parameters in
   `workflow/README.md`). The 610-haplotype run was staged: extraction (results in `mhc_fastas/`), Immuannot
   (`immuannot_gtfs/`), then `nig/stage3.sh mhc_fastas immuannot_gtfs results-stage4` for aggregation, gene cuts,
   bundles and graphs (`results/`).
4. Class II span: `workflow/scripts/extract_span.py --from-gene HLA-DRA --to-gene HLA-DMA`, then
   `toil-cwl-runner workflow/tools/pgr_pbundle.cwl classII/job.yml`.
5. Read support: `read_support/fetch.sh` (1000G CRAM regions), `read_support/run.sh` (realignment, pileups),
   `workflow/scripts/novel_alleles.py` (Immuannot `cds_mut` to contig coordinates).
6. Figures: `workflow/analysis/*.py`, run in a directory holding `hla_calls.tsv`, `gene_copy_number.tsv`, `data/` (tables above).

## Caveats

- Novel alleles: only substitution-type differences with read or independent-assembly support should be trusted; indel-only
  "new" alleles are typical HiFi homopolymer errors. APR and JaSaPaGe raw reads are not public, so their novel alleles are
  assembly-only.
- pgr-query homology fetch of genes (alternative `gene_source: pgr-query`) needs w=24/k=32 (class I) or w=16/k=24 (DRB1).
- HPRC-Jewish is one individual (HG002); treat its frequencies as anecdotal.
- K-PanRef and CPC haplotypes are graph paths: Minigraph-Cactus clips unaligned sequence, so large insertions private to a
  haplotype can be shortened compared with the original assembly.
