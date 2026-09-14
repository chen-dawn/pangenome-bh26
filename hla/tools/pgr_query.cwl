#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: Fetch a region of interest from a set of PanSN-named sequences with pgr-query
doc: |
  Uses the concatenated sequence set (e.g. all MHC haplotypes) as an in-memory
  pgr-tk database and pulls out the segments homologous to the query (e.g. one
  HLA gene +/- flank). Output records are renamed sample#haplotype#<query name>.
requirements:
  InitialWorkDirRequirement:
    listing:
      - entryname: pgr_query_rename.py
        entry: {$include: ../scripts/pgr_query_rename.py}
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
  EnvVarRequirement:
    envDef:
      RAYON_NUM_THREADS: $(String(inputs.threads))
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/pgr-tk:0.5.1--py38hfa1e82d_1
  SoftwareRequirement:
    packages:
      pgr-tk: {version: ["0.5.1"], specs: ["https://github.com/GeneDx/pgr-tk"]}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 12000
inputs:
  database: {type: File, doc: "PanSN-named fasta used as the pgr-tk sequence database"}
  query: {type: File, doc: "single-record fasta; the record name (nameroot of the file) labels the output"}
  threads: {type: int, default: 4}
  w: {type: int, default: 48}
  k: {type: int, default: 56}
  r: {type: int, default: 2}
  min_span: {type: int, default: 12}
  merge_range_tol: {type: int, default: 20000}
  min_piece: {type: int, default: 500}
  min_anchors: {type: int, default: 5}
  min_query_frac: {type: float, default: 0.5}
baseCommand: [bash, -c]
arguments:
  - valueFrom: |
      set -euo pipefail
      G=$(inputs.query.nameroot)
      pgr-query --fastx-file -w $(inputs.w) -k $(inputs.k) -r $(inputs.r) --min-span $(inputs.min_span) \
        --merge-range-tol $(inputs.merge_range_tol) "$(inputs.database.path)" "$(inputs.query.path)" "$G"_pgrq
      python3 pgr_query_rename.py --hit "$G"_pgrq.000.hit --fasta "$G"_pgrq.000.fa --gene "$G" \
        --min-piece $(inputs.min_piece) --min-anchors $(inputs.min_anchors) --min-query-frac $(inputs.min_query_frac) \
        --out-fasta $G.fa --out-tsv $G.hits.tsv
outputs:
  fasta: {type: File, outputBinding: {glob: "$(inputs.query.nameroot).fa"}}
  hits: {type: File, outputBinding: {glob: "$(inputs.query.nameroot).hits.tsv"}}
