#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: odgi layout + draw 2D rendering of a pangenome graph
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
    coresMin: $(inputs.threads)
    ramMin: 8000
inputs:
  graph: {type: File, doc: "odgi .og"}
  threads: {type: int, default: 4}
baseCommand: [bash, -c]
arguments:
  - valueFrom: |
      set -euo pipefail
      odgi layout -i "$(inputs.graph.path)" -o $(inputs.graph.nameroot).lay -T $(inputs.graph.nameroot).lay.tsv -t $(inputs.threads) -P
      odgi draw -i "$(inputs.graph.path)" -c $(inputs.graph.nameroot).lay -p $(inputs.graph.nameroot).draw.png -C -w 20 -H 1500
outputs:
  layout:
    type: File
    outputBinding: {glob: "$(inputs.graph.nameroot).lay"}
  draw:
    type: File
    outputBinding: {glob: "$(inputs.graph.nameroot).draw.png"}
