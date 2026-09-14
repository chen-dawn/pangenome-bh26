#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool
label: Extract the MHC region from one haplotype assembly with pgr-tk (pgr-query)
doc: |
  pgr-query builds a SHIMMER index of the whole assembly in memory (~9 GB, ~30 s
  for a 3 Gb haplotype with 8 threads) and fetches the target segments homologous
  to the reference MHC query sequence. Hits are filtered (span, anchor count) and
  renamed to PanSN sample#haplotype#contig:start-end; minus-strand hits are
  already reverse-complemented by pgr-query (suffix _rc).
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
      python: {version: ["3"]}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 14000
inputs:
  assembly: {type: File}
  query: {type: File, doc: "fasta with the reference MHC region (+flank)"}
  sample: {type: string}
  haplotype: {type: string}
  threads: {type: int, default: 8}
  w: {type: int, default: 80}
  k: {type: int, default: 56}
  r: {type: int, default: 4}
  min_span: {type: int, default: 64}
  merge_range_tol: {type: int, default: 100000}
  min_piece: {type: int, default: 20000}
  min_anchors: {type: int, default: 20}
baseCommand: [bash, -c]
arguments:
  - valueFrom: |
      set -euo pipefail
      P=$(inputs.sample)_$(inputs.haplotype)
      pgr-query --fastx-file -w $(inputs.w) -k $(inputs.k) -r $(inputs.r) --min-span $(inputs.min_span) \
        --merge-range-tol $(inputs.merge_range_tol) "$(inputs.assembly.path)" "$(inputs.query.path)" "$P"_pgrq
      python3 pgr_query_rename.py --hit "$P"_pgrq.000.hit --fasta "$P"_pgrq.000.fa \
        --sample "$(inputs.sample)" --haplotype "$(inputs.haplotype)" \
        --min-piece $(inputs.min_piece) --min-anchors $(inputs.min_anchors) \
        --out-fasta $P.mhc.fa --out-tsv $P.mhc.tsv
      mv "$P"_pgrq.000.hit "$P"_pgrq.hit
outputs:
  mhc_fasta: {type: File, outputBinding: {glob: "*.mhc.fa"}}
  mhc_tsv: {type: File, outputBinding: {glob: "*.mhc.tsv"}}
  hits: {type: File, outputBinding: {glob: "*_pgrq.hit"}}
