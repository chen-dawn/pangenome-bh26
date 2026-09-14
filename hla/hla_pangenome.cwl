#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: Workflow
label: HLA region extraction, allele annotation and pangenome visualisation
doc: |
  For every haplotype assembly: fetch the MHC region with pgr-tk (pgr-query
  against the GRCh38 MHC; minimap2 projection available as alternative),
  annotate HLA/KIR/C4 genes and call alleles (Immuannot, IPD-IMGT/HLA), then
  aggregate calls across cohorts, fetch every HLA gene from all MHC haplotypes
  with pgr-query, and visualise each gene as a pgr-tk principal-bundle plot and
  as a pggb graph rendered with odgi viz / odgi draw. The whole MHC set gets a
  pgr-tk bundle decomposition too; a whole-MHC pggb graph is optional.
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  SubworkflowFeatureRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  assemblies: {type: "File[]", doc: "haplotype-resolved assemblies (fa or fa.gz)"}
  samples: {type: "string[]", doc: "sample id per assembly"}
  haplotypes: {type: "string[]", doc: "haplotype label per assembly (1/2)"}
  cohorts: {type: "string[]", doc: "cohort label per assembly"}
  reference: {type: File, doc: "GRCh38 fasta with PanSN contig names"}
  reference_chrom: {type: string, default: "GRCh38#0#chr6"}
  mhc_region: {type: string, default: "GRCh38#0#chr6:28510120-33480577", doc: "extended MHC on GRCh38"}
  mhc_flank: {type: int, default: 100000}
  extractor: {type: string, default: "pgr-tk", doc: "pgr-tk (pgr-query) or minimap2"}
  immuannot_dir: {type: Directory}
  immuannot_ref: {type: Directory}
  graph_genes: {type: string, default: "HLA-A,HLA-B,HLA-C,HLA-E,HLA-F,HLA-G,HLA-DRA,HLA-DRB1,HLA-DRB3,HLA-DRB4,HLA-DRB5,HLA-DQA1,HLA-DQB1,HLA-DPA1,HLA-DPB1,MICA,MICB,TAP1,TAP2,C4A,C4B"}
  gene_flank: {type: int, default: 2000}
  extract_threads: {type: int, default: 8}
  immuannot_threads: {type: int, default: 4}
  pggb_threads: {type: int, default: 8}
  gene_percent_identity: {type: int, default: 90}
  gene_segment_length: {type: int, default: 2000}
  gene_bundle_w: {type: int, default: 48}
  gene_bundle_r: {type: int, default: 2}
  gene_bundle_length_cutoff: {type: int, default: 200}
  gene_bundle_merge_distance: {type: int, default: 1000}
  mhc_bundle_w: {type: int, default: 80}
  mhc_bundle_r: {type: int, default: 6}
  mhc_bundle_min_span: {type: int, default: 28}
  build_mhc_graph: {type: boolean, default: false, doc: "also build a whole-MHC pggb graph (slow for hundreds of haplotypes)"}
  mhc_graph_n: {type: int, default: 24, doc: "pggb -n for the whole-MHC graph (sparsifies mapping)"}
  mhc_graph_threads: {type: int, default: 24}

steps:
  prepare_reference:
    run: tools/prepare_reference.cwl
    in:
      reference: reference
      chrom: reference_chrom
      region: mhc_region
      flank: mhc_flank
    out: [chrom_fasta, chrom_index, mhc_query]

  extract_pgrtk:
    run: tools/extract_mhc_pgrtk.cwl
    when: $(inputs.extractor == "pgr-tk")
    scatter: [assembly, sample, haplotype]
    scatterMethod: dotproduct
    in:
      extractor: extractor
      assembly: assemblies
      sample: samples
      haplotype: haplotypes
      query: prepare_reference/mhc_query
      threads: extract_threads
    out: [mhc_fasta, mhc_tsv]

  extract_minimap2:
    run: tools/extract_mhc.cwl
    when: $(inputs.extractor == "minimap2")
    scatter: [assembly, sample, haplotype]
    scatterMethod: dotproduct
    in:
      extractor: extractor
      assembly: assemblies
      sample: samples
      haplotype: haplotypes
      reference_index: prepare_reference/chrom_index
      region: mhc_region
      flank: mhc_flank
      threads: extract_threads
    out: [mhc_fasta, mhc_tsv]

  immuannot:
    run: tools/immuannot.cwl
    scatter: [contigs, sample, haplotype]
    scatterMethod: dotproduct
    in:
      contigs:
        source: [extract_pgrtk/mhc_fasta, extract_minimap2/mhc_fasta]
        pickValue: first_non_null
      sample: samples
      haplotype: haplotypes
      label:
        valueFrom: $(inputs.sample)_$(inputs.haplotype)
      immuannot_dir: immuannot_dir
      immuannot_ref: immuannot_ref
      threads: immuannot_threads
    out: [gtf, log]

  aggregate:
    run: tools/aggregate_immuannot.cwl
    in:
      gtfs: immuannot/gtf
      mhc_fastas:
        source: [extract_pgrtk/mhc_fasta, extract_minimap2/mhc_fasta]
        pickValue: first_non_null
      mhc_tsvs:
        source: [extract_pgrtk/mhc_tsv, extract_minimap2/mhc_tsv]
        pickValue: first_non_null
      samples: samples
      haplotypes: haplotypes
      cohorts: cohorts
      genes: graph_genes
      gene_flank: gene_flank
    out: [calls, calls_matrix, copy_number, extraction_summary, gene_queries, gene_list, plots]

  concat_mhc:
    run:
      class: CommandLineTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        fastas: {type: "File[]"}
      baseCommand: cat
      arguments: [$(inputs.fastas)]
      stdout: MHC.fa
      outputs:
        merged: {type: File, outputBinding: {glob: MHC.fa}}
    in:
      fastas:
        source: [extract_pgrtk/mhc_fasta, extract_minimap2/mhc_fasta]
        pickValue: first_non_null
    out: [merged]

  gene_fetch:
    run: tools/pgr_query.cwl
    scatter: query
    in:
      database: concat_mhc/merged
      query: aggregate/gene_queries
    out: [fasta, hits]

  gene_bundle:
    run: tools/pgr_pbundle.cwl
    scatter: sequences
    in:
      sequences: gene_fetch/fasta
      w: gene_bundle_w
      r: gene_bundle_r
      bundle_length_cutoff: gene_bundle_length_cutoff
      bundle_merge_distance: gene_bundle_merge_distance
    out: [bed, summary, mapg_gfa, pmapg_gfa, dist, newick, order, svg, html]

  gene_graph:
    run: tools/pggb.cwl
    scatter: sequences
    in:
      sequences: gene_fetch/fasta
      percent_identity: gene_percent_identity
      segment_length: gene_segment_length
      threads: pggb_threads
    out: [gfa, og, stats, pggb_dir]

  gene_viz:
    run: tools/odgi_viz.cwl
    scatter: graph
    in:
      graph: gene_graph/og
    out: [viz, viz_depth]

  gene_draw:
    run: tools/odgi_draw.cwl
    scatter: graph
    in:
      graph: gene_graph/og
    out: [draw, layout]

  mhc_bundle:
    run: tools/pgr_pbundle.cwl
    in:
      sequences: concat_mhc/merged
      w: mhc_bundle_w
      r: mhc_bundle_r
      min_span: mhc_bundle_min_span
      threads: {default: 8}
    out: [bed, summary, mapg_gfa, pmapg_gfa, dist, newick, order, svg, html]

  mhc_graph:
    run: tools/pggb.cwl
    when: $(inputs.enabled)
    in:
      enabled: build_mhc_graph
      sequences: concat_mhc/merged
      n_haplotypes: mhc_graph_n
      segment_length: {default: 10000}
      threads: mhc_graph_threads
    out: [gfa, og, stats, pggb_dir]

  mhc_viz:
    run: tools/odgi_viz.cwl
    when: $(inputs.enabled)
    in:
      enabled: build_mhc_graph
      graph: mhc_graph/og
      height: {default: 1500}
    out: [viz, viz_depth]

outputs:
  mhc_fastas:
    type: "File[]"
    outputSource: [extract_pgrtk/mhc_fasta, extract_minimap2/mhc_fasta]
    pickValue: first_non_null
  mhc_all: {type: File, outputSource: concat_mhc/merged}
  immuannot_gtfs: {type: "File[]", outputSource: immuannot/gtf}
  immuannot_logs: {type: "File[]", outputSource: immuannot/log}
  hla_calls: {type: File, outputSource: aggregate/calls}
  hla_calls_matrix: {type: File, outputSource: aggregate/calls_matrix}
  gene_copy_number: {type: File, outputSource: aggregate/copy_number}
  mhc_extraction_summary: {type: File, outputSource: aggregate/extraction_summary}
  gene_list: {type: File, outputSource: aggregate/gene_list}
  gene_queries: {type: "File[]", outputSource: aggregate/gene_queries}
  summary_plots: {type: "File[]", outputSource: aggregate/plots}
  gene_fastas: {type: "File[]", outputSource: gene_fetch/fasta}
  gene_hits: {type: "File[]", outputSource: gene_fetch/hits}
  gene_bundle_svg: {type: "File[]", outputSource: gene_bundle/svg}
  gene_bundle_html: {type: "File[]", outputSource: gene_bundle/html}
  gene_bundle_bed: {type: "File[]", outputSource: gene_bundle/bed}
  gene_bundle_gfa: {type: "File[]", outputSource: gene_bundle/pmapg_gfa}
  gene_bundle_newick: {type: "File[]", outputSource: gene_bundle/newick}
  gene_bundle_summary: {type: "File[]", outputSource: gene_bundle/summary}
  gene_gfas: {type: "File[]", outputSource: gene_graph/gfa}
  gene_ogs: {type: "File[]", outputSource: gene_graph/og}
  gene_graph_stats: {type: "File[]", outputSource: gene_graph/stats}
  gene_viz_png: {type: "File[]", outputSource: gene_viz/viz}
  gene_viz_depth_png: {type: "File[]?", outputSource: gene_viz/viz_depth}
  gene_draw_png: {type: "File[]", outputSource: gene_draw/draw}
  mhc_bundle_svg: {type: File, outputSource: mhc_bundle/svg}
  mhc_bundle_html: {type: File, outputSource: mhc_bundle/html}
  mhc_bundle_bed: {type: File, outputSource: mhc_bundle/bed}
  mhc_bundle_gfa: {type: File, outputSource: mhc_bundle/pmapg_gfa}
  mhc_bundle_newick: {type: File, outputSource: mhc_bundle/newick}
  mhc_bundle_summary: {type: File, outputSource: mhc_bundle/summary}
  mhc_gfa: {type: File?, outputSource: mhc_graph/gfa}
  mhc_og: {type: File?, outputSource: mhc_graph/og}
  mhc_graph_stats: {type: File?, outputSource: mhc_graph/stats}
  mhc_viz_png: {type: File?, outputSource: mhc_viz/viz}
