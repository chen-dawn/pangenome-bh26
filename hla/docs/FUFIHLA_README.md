# FuFiHLA typing of the assembled MHC haplotypes (`/home/asianhla/data/upload/HLA/fufihla/`)

FuFiHLA 0.2.4 (Hu et al., "FuFiHLA: a tool for full-field HLA typing from long-read data", Bioinformatics 42(5) btag231,
2026; https://github.com/jingqing-hu/FuFiHLA) types HLA-A, -B, -C, -DRB1, -DQA1 and -DQB1 at full field resolution from
PacBio HiFi or ONT reads. The assemblies used here have no public raw long reads (APR, JaSaPaGe, K-PanRef, CPC), so FuFiHLA
is run on the assembled sequence itself:

- both assembled MHC haplotypes of an individual (`<sample>_1.mhc.fa`, `<sample>_2.mhc.fa`, the same MHC fastas that
  Immuannot annotated) are tiled into error-free 15 kb pseudo-HiFi reads every 1.5 kb (10x per haplotype;
  `workflow/scripts/tile_reads.py`), and the combined diploid read set is typed with
  `fufihla --fa reads.fa.gz --out <sample> --refdir ref_data` (reference alleles from `fufihla-ref-prep`:
  IPD-IMGT/HLA 3.65.0, whereas Immuannot's database is IPD-IMGT/HLA 3.55; many sequences that Immuannot reports as
  `:new` have since been named, e.g. DRB1*09:59, DRB1*14:06:01, A*24:608N)
- check on HG002: all 12 alleles identical to Immuannot at full resolution (including the novel DRB1*04:02:01 variant)
- reference genomes (GRCh38, CHM13; single haplotype) are not typed

## Files

- `individuals.tsv` - sample, hap1 MHC fasta, hap2 MHC fasta (array index = line number)
- `calls/<sample>.out` - FuFiHLA final calls (column 1 allele, column 2 consensus id `cons_h1/h2_<allele>`, last column
  minimap2 `cs` string of the consensus vs the allele: `cs:Z::<len>` = exact match, otherwise novel differences)
- `calls/<sample>.log`, `consensus/<sample>/` - run log and reconstructed allele consensus sequences
- `logs/` - Slurm array logs
- `../tables/typing_concordance*.tsv` - Immuannot vs FuFiHLA vs T1K (1000G Illumina reads, `../../1000G_MHC/`)

## How it was generated

```bash
micromamba create -n typing -c conda-forge -c bioconda fufihla t1k samtools
fufihla-ref-prep                                   # in ~/hla/typing -> ref_data/
sbatch -p epyc,medium -A general_analysis --array=1-N%40 workflow/nig/fufihla_assemblies.sbatch
```

## Caveats

- Pseudo-reads from an assembly carry the assembly's errors; FuFiHLA is then a second, independent allele-calling
  algorithm on the same sequence, not an independent measurement. Disagreements with Immuannot point at calling
  differences (template choice, novel alleles), not at sequencing error.
- An individual whose two assembled haplotypes are identical (JaSaPaGe NA18952) is typed homozygous by construction.
