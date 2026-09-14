#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: Extract the MHC region from one haplotype assembly
doc: |
  minimap2 (asm20) maps the assembly contigs to the reference chromosome; the
  reference MHC interval (+flank) is projected through the alignments onto the
  contigs and the covering segments are written in reference orientation with
  PanSN names sample#haplotype#contig:start-end.
requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: extract_mhc.py
        entry: {$include: ../scripts/extract_mhc.py}
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/minimap2:2.28--he4a0461_3
  SoftwareRequirement:
    packages:
      minimap2: {version: ["2.28"]}
      python: {version: ["3"]}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 12000
inputs:
  assembly: {type: File}
  reference_index: {type: File, doc: "minimap2 .mmi of the single reference chromosome"}
  region: {type: string, doc: "chrom:start-end of the MHC on the reference chromosome"}
  sample: {type: string}
  haplotype: {type: string}
  flank: {type: int, default: 100000}
  threads: {type: int, default: 4}
  preset: {type: string, default: asm20}
baseCommand: [python3, extract_mhc.py]
arguments:
  - {prefix: --assembly, valueFrom: $(inputs.assembly.path)}
  - {prefix: --reference, valueFrom: $(inputs.reference_index.path)}
  - {prefix: --region, valueFrom: $(inputs.region)}
  - {prefix: --sample, valueFrom: $(inputs.sample)}
  - {prefix: --haplotype, valueFrom: $(inputs.haplotype)}
  - {prefix: --flank, valueFrom: $(inputs.flank)}
  - {prefix: --threads, valueFrom: $(inputs.threads)}
  - {prefix: --preset, valueFrom: $(inputs.preset)}
  - {prefix: --out-fasta, valueFrom: $(inputs.sample)_$(inputs.haplotype).mhc.fa}
  - {prefix: --out-tsv, valueFrom: $(inputs.sample)_$(inputs.haplotype).mhc.tsv}
  - {prefix: --out-paf, valueFrom: $(inputs.sample)_$(inputs.haplotype).chr6.paf}
outputs:
  mhc_fasta:
    type: File
    outputBinding: {glob: "*.mhc.fa"}
  mhc_tsv:
    type: File
    outputBinding: {glob: "*.mhc.tsv"}
  paf:
    type: File
    outputBinding: {glob: "*.chr6.paf"}
