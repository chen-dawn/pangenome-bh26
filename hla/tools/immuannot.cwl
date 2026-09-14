#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: Immuannot HLA/KIR/C4 gene annotation and allele calling on assembly contigs
doc: |
  Immuannot (Zhou et al., Genome Research 2024, doi:10.1101/gr.278985.124)
  aligns IPD-IMGT/HLA, IPD-KIR and RefSeq gene sequences to the contigs with
  minimap2 and reports gene structure plus full-resolution allele calls in a
  gzip GTF. The tool code and the reference bundle (Zenodo 10.5281/zenodo.10948964,
  Data-2024Feb02) are passed in as directories so versions are explicit.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/minimap2:2.28--he4a0461_3
  SoftwareRequirement:
    packages:
      minimap2: {version: ["2.28"]}
      python: {version: ["3"]}
      immuannot: {specs: ["https://github.com/YingZhou001/Immuannot"]}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 6000
inputs:
  contigs: {type: File}
  immuannot_dir: {type: Directory, doc: "clone of the Immuannot repository"}
  immuannot_ref: {type: Directory, doc: "extracted Data-YYYYMMMDD reference bundle"}
  label: {type: string}
  threads: {type: int, default: 4}
  diff: {type: float, default: 0.03}
baseCommand: [bash]
arguments:
  - valueFrom: $(inputs.immuannot_dir.path)/scripts.pub.v3/immuannot.sh
  - {prefix: -c, valueFrom: $(inputs.contigs.path)}
  - {prefix: -r, valueFrom: $(inputs.immuannot_ref.path)}
  - {prefix: -o, valueFrom: $(inputs.label)}
  - {prefix: -t, valueFrom: $(inputs.threads)}
  - {prefix: --diff, valueFrom: $(inputs.diff)}
stdout: $(inputs.label).immuannot.log
outputs:
  gtf:
    type: File
    outputBinding: {glob: "$(inputs.label).gtf.gz"}
  log:
    type: File
    outputBinding: {glob: "$(inputs.label).immuannot.log"}
