#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: Cut annotated HLA genes out of the MHC haplotypes (Immuannot coordinates)
doc: |
  For every gene in the list, extracts each annotated copy (gene +/- flank, in
  gene orientation) from the concatenated MHC fasta with samtools faidx, using
  the coordinates Immuannot reported (hla_calls.tsv). Records are named
  sample#haplotype#GENE[_n]. This recovers every annotated haplotype, unlike a
  homology search seeded from a single reference allele.
requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: extract_genes.py
        entry: {$include: ../scripts/extract_genes.py}
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/samtools:1.21--h50ea8bc_0
  SoftwareRequirement:
    packages:
      samtools: {version: ["1.21"]}
      python: {version: ["3"]}
  ResourceRequirement:
    coresMin: 1
    ramMin: 4000
inputs:
  calls: {type: File, doc: "hla_calls.tsv from aggregate_immuannot"}
  fasta: {type: File, doc: "concatenated MHC haplotypes (PanSN names)"}
  genes: {type: string}
  flank: {type: int, default: 2000}
  min_seqs: {type: int, default: 4}
baseCommand: [python3, extract_genes.py]
arguments:
  - {prefix: --calls, valueFrom: $(inputs.calls)}
  - {prefix: --fasta, valueFrom: $(inputs.fasta)}
  - {prefix: --genes, valueFrom: $(inputs.genes)}
  - {prefix: --flank, valueFrom: $(inputs.flank)}
  - {prefix: --min-seqs, valueFrom: $(inputs.min_seqs)}
  - {prefix: --outdir, valueFrom: genes}
outputs:
  fastas: {type: "File[]", outputBinding: {glob: "genes/*.fa"}}
  regions: {type: "File[]", outputBinding: {glob: "genes/*.regions.tsv"}}
