#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: Aggregate Immuannot calls, write per-gene query sequences, plot allele frequencies
requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: aggregate_immuannot.py
        entry: {$include: ../scripts/aggregate_immuannot.py}
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/mulled-v2-f42a44964bca5225c7860882c231a7db2f8dc5ea:2d6a5ac91cbc2a5ee4c9d9e6cb6d0f8c1d47b7bf-0
  SoftwareRequirement:
    packages:
      python: {version: ["3.11"]}
      pandas: {}
      matplotlib: {}
  ResourceRequirement:
    coresMin: 2
    ramMin: 16000
inputs:
  gtfs: {type: "File[]"}
  mhc_fastas: {type: "File[]"}
  mhc_tsvs: {type: "File[]"}
  samples: {type: "string[]"}
  haplotypes: {type: "string[]"}
  cohorts: {type: "string[]"}
  genes: {type: string, default: "HLA-A,HLA-B,HLA-C,HLA-E,HLA-F,HLA-G,HLA-DRA,HLA-DRB1,HLA-DRB3,HLA-DRB4,HLA-DRB5,HLA-DQA1,HLA-DQB1,HLA-DPA1,HLA-DPB1,MICA,MICB,TAP1,TAP2,C4A,C4B"}
  gene_flank: {type: int, default: 2000}
  min_graph_seqs: {type: int, default: 4}
baseCommand: [python3, aggregate_immuannot.py]
arguments:
  - {prefix: --gtf, valueFrom: $(inputs.gtfs)}
  - {prefix: --mhc-fasta, valueFrom: $(inputs.mhc_fastas)}
  - {prefix: --mhc-tsv, valueFrom: $(inputs.mhc_tsvs)}
  - {prefix: --sample, valueFrom: $(inputs.samples)}
  - {prefix: --haplotype, valueFrom: $(inputs.haplotypes)}
  - {prefix: --cohort, valueFrom: $(inputs.cohorts)}
  - {prefix: --genes, valueFrom: $(inputs.genes)}
  - {prefix: --gene-flank, valueFrom: $(inputs.gene_flank)}
  - {prefix: --min-graph-seqs, valueFrom: $(inputs.min_graph_seqs)}
  - {prefix: --outdir, valueFrom: "."}
outputs:
  calls: {type: File, outputBinding: {glob: hla_calls.tsv}}
  calls_matrix: {type: File, outputBinding: {glob: hla_calls_matrix.tsv}}
  copy_number: {type: File, outputBinding: {glob: gene_copy_number.tsv}}
  extraction_summary: {type: File, outputBinding: {glob: mhc_extraction_summary.tsv}}
  gene_queries: {type: "File[]", outputBinding: {glob: "queries/*.fa"}}
  gene_list: {type: File, outputBinding: {glob: genes.txt}}
  plots: {type: "File[]", outputBinding: {glob: "plots/*.png"}}
