#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: Build a pangenome graph with pggb
doc: |
  PanGenome Graph Builder (wfmash + seqwish + smoothxg + gfaffix + odgi).
  Input fasta must use PanSN names; it is bgzipped and indexed here.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: ghcr.io/pangenome/pggb:latest
  SoftwareRequirement:
    packages:
      pggb: {version: ["0.7.4"]}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 16000
inputs:
  sequences: {type: File, doc: "PanSN-named fasta; output files are named after its nameroot"}
  n_haplotypes: {type: int?, doc: "pggb -n (mappings per segment); defaults to the number of sequences"}
  percent_identity: {type: int, default: 90}
  segment_length: {type: int, default: 5000}
  poa_params: {type: string?, default: "asm20"}
  threads: {type: int, default: 8}
baseCommand: [bash, -c]
arguments:
  - valueFrom: |
      set -euo pipefail
      NAME=$(inputs.sequences.nameroot)
      N=$(inputs.n_haplotypes || 0)
      cp "$(inputs.sequences.path)" $NAME.fa
      [ "$N" -gt 0 ] || N=`grep -c '>' $NAME.fa`
      bgzip -@ $(inputs.threads) $NAME.fa
      samtools faidx $NAME.fa.gz
      pggb -i $NAME.fa.gz -o pggb_$NAME \
        -n $N -p $(inputs.percent_identity) -s $(inputs.segment_length) \
        -P $(inputs.poa_params) -t $(inputs.threads) -S 2>&1 | tail -50
      cp pggb_$NAME/*.smooth.final.gfa $NAME.gfa
      cp pggb_$NAME/*.smooth.final.og $NAME.og
      odgi stats -i $NAME.og -S > $NAME.stats.tsv
outputs:
  gfa:
    type: File
    outputBinding: {glob: "$(inputs.sequences.nameroot).gfa"}
  og:
    type: File
    outputBinding: {glob: "$(inputs.sequences.nameroot).og"}
  stats:
    type: File
    outputBinding: {glob: "$(inputs.sequences.nameroot).stats.tsv"}
  pggb_dir:
    type: Directory
    outputBinding: {glob: "pggb_$(inputs.sequences.nameroot)"}
