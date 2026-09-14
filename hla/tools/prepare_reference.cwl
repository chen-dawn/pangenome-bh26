#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: Extract one reference chromosome and build a minimap2 index for it
doc: |
  Pulls a single chromosome (default GRCh38 chr6, PanSN name GRCh38#0#chr6)
  out of the full reference and builds a minimap2 asm-preset index so every
  scattered extraction step can map against the small target quickly.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0
  SoftwareRequirement:
    packages:
      samtools: {version: ["1.21"]}
      minimap2: {version: ["2.28"]}
  ResourceRequirement:
    coresMin: 2
    ramMin: 8000
inputs:
  reference:
    type: File
    doc: whole-genome reference fasta (uncompressed or bgzip)
  chrom:
    type: string
    default: GRCh38#0#chr6
  preset:
    type: string
    default: asm20
  region:
    type: string
    default: "GRCh38#0#chr6:28510120-33480577"
    doc: "MHC interval on the reference chromosome (1-based, inclusive)"
  flank:
    type: int
    default: 100000
baseCommand: [bash, -c]
arguments:
  - valueFrom: |
      set -euo pipefail
      samtools faidx "$(inputs.reference.path)" "$(inputs.chrom)" > chrom.fa
      samtools faidx chrom.fa
      minimap2 -x $(inputs.preset) -d chrom.mmi chrom.fa
      CHROM=`echo "$(inputs.region)" | sed 's/:[^:]*$//'`
      S=`echo "$(inputs.region)" | sed 's/.*://; s/-.*//'`
      E=`echo "$(inputs.region)" | sed 's/.*-//'`
      S=`expr $S - $(inputs.flank)`; [ "$S" -ge 1 ] || S=1
      E=`expr $E + $(inputs.flank)`
      samtools faidx chrom.fa "$CHROM:$S-$E" | sed "1s/.*/>MHC_$CHROM:$S-$E/" > mhc_query.fa
outputs:
  chrom_fasta:
    type: File
    outputBinding: {glob: chrom.fa}
    secondaryFiles: [.fai]
  chrom_index:
    type: File
    outputBinding: {glob: chrom.mmi}
  mhc_query:
    type: File
    outputBinding: {glob: mhc_query.fa}
