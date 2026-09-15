#!/bin/bash
# Regenerate every table and figure on the cluster after a workflow run (no local computation).
# usage: bash nig/make_figures.sh <results dir, e.g. ~/hla/results-stage5> <classII dir, e.g. ~/hla/classII5> <viz dir, e.g. ~/hla/viz5>
# viz dir must already hold data/ with the read-support inputs (rs/, novel_coding_support.tsv, het_density_100kb.tsv,
# mhc_9samples.gt.tsv) from the earlier read-support analysis; everything else is (re)built here.
set -euo pipefail
R=${1:?results}; C=${2:?classII}; V=${3:?viz}
CWL=$HOME/hla/cwl
export PATH=$HOME/hla/mm/envs/hla/bin:$HOME/hla/mm/envs/typing/bin:$PATH OPENBLAS_NUM_THREADS=1 PYTHONPATH=$CWL/analysis
mkdir -p $V/data/classII $V/figures $V/gene_bundles $V/gene_graphs $V/plots $V/mhc_bundle $V/classII_bundle
cd $V
cp $R/hla_calls.tsv $R/gene_copy_number.tsv $R/hla_calls_matrix.tsv $R/mhc_extraction_summary.tsv .
cp $CWL/data/hprc_r2_populations.tsv data/
cp $C/out/classII.* data/classII/ && cp $C/out/classII.html $C/out/classII.svg classII_bundle/
cp $R/MHC.html $R/MHC.svg mhc_bundle/ 2>/dev/null || true
for g in $(cat $R/genes.txt); do cp $R/$g.html $R/$g.svg gene_bundles/ 2>/dev/null || true; cp $R/$g.viz.png $R/$g.viz_depth.png $R/$g.draw.png gene_graphs/ 2>/dev/null || true; done
cp $R/*.png plots/ 2>/dev/null || true
# novel coding differences from every GTF (all cohorts)
python3 $CWL/scripts/novel_alleles.py $R/*.gtf.gz > data/novel_coding_all.tsv
echo "novel coding rows: $(($(wc -l < data/novel_coding_all.tsv)-1))"
python3 $CWL/analysis/fig2_extraction.py
python3 $CWL/analysis/fig5_population.py
python3 $CWL/analysis/classII_flow_pca.py
python3 $CWL/analysis/classII_bundles.py
python3 $CWL/analysis/fig3_mhc_homozygosity.py
python3 $CWL/analysis/fig4_novel_alleles.py
python3 $CWL/analysis/typing_concordance.py --calls hla_calls.tsv \
  --fufihla /home/asianhla/data/upload/HLA/fufihla/calls --gene-fasta-dir $R \
  --t1k /home/asianhla/data/upload/1000G_MHC/t1k \
  --published $CWL/ground_truth/1000G_2014/20140702_hla_diversity.txt \
  --hprc-truth $CWL/ground_truth/HPRC_2024/Lai2024_HPRC_Supplementary-6_HLA_4field_genotypes.xlsx \
  --out-prefix data/typing_concordance --fig figures/fig12_typing_concordance.png
ls -la figures
echo FIGSDONE
