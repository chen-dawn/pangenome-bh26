#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: pgr-tk principal bundle decomposition and SVG rendering
doc: |
  pgr-pbundle-decomp builds the MAP-graph of the sequence set and decomposes each
  sequence into principal bundles (bed + GFA); pgr-pbundle-bed2dist derives a
  bundle-based distance matrix and dendrogram; pgr-pbundle-bed2svg renders the
  bundle tracks (one row per haplotype, ordered by the dendrogram) as SVG/HTML.
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  EnvVarRequirement:
    envDef:
      RAYON_NUM_THREADS: $(inputs.threads)
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/pgr-tk:0.5.1--py38hfa1e82d_1
  SoftwareRequirement:
    packages:
      pgr-tk: {version: ["0.5.1"], specs: ["https://github.com/GeneDx/pgr-tk"]}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 16000
inputs:
  sequences: {type: File, doc: "PanSN-named fasta; outputs are named after its nameroot"}
  threads: {type: int, default: 4}
  w: {type: int, default: 48}
  k: {type: int, default: 56}
  r: {type: int, default: 4}
  min_span: {type: int, default: 12}
  min_cov: {type: int, default: 0}
  min_branch_size: {type: int, default: 8}
  bundle_length_cutoff: {type: int, default: 2500}
  bundle_merge_distance: {type: int, default: 10000}
  track_panel_width: {type: int, default: 1600}
  track_tick_interval: {type: int?}
baseCommand: [bash, -c]
arguments:
  - valueFrom: |
      set -euo pipefail
      P=$(inputs.sequences.nameroot)
      pgr-pbundle-decomp -w $(inputs.w) -k $(inputs.k) -r $(inputs.r) --min-span $(inputs.min_span) \
        --min-cov $(inputs.min_cov) --min-branch-size $(inputs.min_branch_size) \
        --bundle-length-cutoff $(inputs.bundle_length_cutoff) --bundle-merge-distance $(inputs.bundle_merge_distance) \
        "$(inputs.sequences.path)" $P
      pgr-pbundle-bed2dist $P.bed $P
      pgr-pbundle-bed2sorted $P.bed $P
      RANGE=`awk 'NR>1 && $2>m {m=$2} END {print m}' $P.ctg.summary.tsv`
      TICK=$(inputs.track_tick_interval || 0)
      [ "$TICK" -gt 0 ] || TICK=`python3 -c "import math; r=$RANGE; e=10**int(math.log10(r/5)); print(int(e*max(1,round(r/5/e))))"`
      pgr-pbundle-bed2svg $P.bed $P --ddg-file $P.ddg --track-range $RANGE --track-tick-interval $TICK \
        --track-panel-width $(inputs.track_panel_width) --html
outputs:
  bed: {type: File, outputBinding: {glob: "$(inputs.sequences.nameroot).bed"}}
  summary: {type: File, outputBinding: {glob: "$(inputs.sequences.nameroot).ctg.summary.tsv"}}
  mapg_gfa: {type: File, outputBinding: {glob: "$(inputs.sequences.nameroot).mapg.gfa"}}
  pmapg_gfa: {type: File, outputBinding: {glob: "$(inputs.sequences.nameroot).pmapg.gfa"}}
  dist: {type: File, outputBinding: {glob: "$(inputs.sequences.nameroot).dist"}}
  newick: {type: File, outputBinding: {glob: "$(inputs.sequences.nameroot).nwk"}}
  order: {type: File, outputBinding: {glob: "$(inputs.sequences.nameroot).ord"}}
  svg: {type: File, outputBinding: {glob: "$(inputs.sequences.nameroot).svg"}}
  html: {type: File, outputBinding: {glob: "$(inputs.sequences.nameroot).html"}}
