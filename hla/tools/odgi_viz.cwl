#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: odgi viz 1D rendering of a pangenome graph
doc: One row per haplotype path, coloured by sample; also a variant-highlighted rendering.
requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/odgi:0.9.4--h077b44d_0
  SoftwareRequirement:
    packages:
      odgi: {version: ["0.9.4"]}
  ResourceRequirement:
    coresMin: 2
    ramMin: 8000
inputs:
  graph: {type: File, doc: "odgi .og"}
  width: {type: int, default: 1800}
  height: {type: int, default: 600}
baseCommand: [bash, -c]
arguments:
  - valueFrom: |
      set -euo pipefail
      odgi viz -i "$(inputs.graph.path)" -o $(inputs.graph.nameroot).viz.png -x $(inputs.width) -y $(inputs.height) -a 10 -s '#' -P
      odgi viz -i "$(inputs.graph.path)" -o $(inputs.graph.nameroot).viz_pos.png -x $(inputs.width) -y $(inputs.height) -a 10 -s '#' -P -du
outputs:
  viz:
    type: File
    outputBinding: {glob: "$(inputs.graph.nameroot).viz.png"}
  viz_pos:
    type: File
    outputBinding: {glob: "$(inputs.graph.nameroot).viz_pos.png"}
