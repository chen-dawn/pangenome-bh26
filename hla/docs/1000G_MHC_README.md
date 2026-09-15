# 1000 Genomes MHC reads and HLA typing (`/home/asianhla/data/upload/1000G_MHC/`)

Generated September 2026 (BioHackathon Japan) by Robert Hoehndorf (KAUST) with Claude Code on the NIG supercomputer.
Scripts: `workflow/nig/kg_mhc_reads.sbatch` (git: https://github.com/leechuck/pangenome-bh26, directory `hla/`).

## Source

- 2504 unrelated 1000 Genomes samples, NYGC 30x Illumina, CRAMs on the NIG shared data area
  `/usr/local/shared_data/public-human-genomes/GRCh38/1000Genomes/CRAM/<sample>/<sample>.cram` (read-only, visible from
  a001 and the epyc nodes), aligned with bwa-mem 0.7.15 to `GRCh38_full_analysis_set_plus_decoy_hla.fa`
  (`/usr/local/shared_data/public-human-genomes/GRCh38/fasta/`; a symlink with its own .fai lives in `~leechuck/hla/ref1kg/`
  because samtools needs a writable place for the index).

## What is here

- `samples.txt` - the 2504 sample IDs (array index = line number)
- `regions.txt` - extracted regions: GRCh38 `chr6:28510120-33480577` (the MHC region used for all assembly extraction),
  all 16 `chr6_*_alt` contigs (includes the MHC alternative haplotypes GL000250-GL000256) and the 525 `HLA-*` allele contigs
- `cram/<sample>.mhc.cram` (+ `.crai`) - every read aligned to those regions, as CRAM against the same reference
  (about 70 MB and 1.7 M reads per sample). Decode with `samtools view -T GRCh38_full_analysis_set_plus_decoy_hla.fa`.
- `t1k/<sample>/` - HLA typing with T1K (Song et al., Genome Research 2023; bioconda t1k) on these reads:
  `<sample>_genotype.tsv` (gene, number of alleles, allele 1, abundance, quality, allele 2, abundance, quality),
  `_allele.tsv`, `_candidate`/aligned read files. Index: IPD-IMGT/HLA downloaded on 2026-09-15
  (`t1k-build.pl -o hlaidx --download IPD-IMGT/HLA`, genomic sequences `hlaidx_dna_seq.fa`).
- `published/20140702_hla_diversity.txt` - published 1000G HLA genotypes (Gourraud et al., PLoS One 2014; Sanger-based
  and PCR-SSO typing of 1000G phase 1-3 cell lines), used as external truth set for T1K
- `logs/` - Slurm array logs

## How it was generated

```bash
sbatch -p epyc -A general_analysis --array=1-2504%40 workflow/nig/kg_mhc_reads.sbatch
```

Per sample: `samtools view -C -T <ref> <cram> $(cat regions.txt)`; `samtools collate | samtools fastq -F 0x900` (primary
reads, pairs); `run-t1k -1 r1.fq.gz -2 r2.fq.gz --preset hla -f hlaidx_dna_seq.fa -t 4`.

## Relation to the assembly-based typing

The same individuals that were assembled (HPRC r2 and the JaSaPaGe Japanese samples are 1000G individuals) are typed three
ways: Immuannot on the assembled haplotypes, FuFiHLA on tiled pseudo-reads of the assembled haplotypes
(`../HLA/fufihla/`), and T1K on these Illumina reads. `workflow/analysis/typing_concordance.py` compares them.

## Caveats

- Reads from HLA-like sequence outside the extracted regions (e.g. unmapped reads, other decoys) are not included.
- T1K resolves alleles from short reads; field depth beyond 2-3 fields is not reliable for all genes.
