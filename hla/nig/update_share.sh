#!/bin/bash
# Refresh the team folder /home/asianhla/data/upload/HLA after the 754-haplotype rerun (stage 5).
# Keeps the previous 610-haplotype results as results-610/ (never deleted), copies new results, class II, MHC/GTF
# folders, figures, tables and slides, and makes everything group-readable.
set -euo pipefail
H=$HOME/hla; S=/home/asianhla/data/upload/HLA
[ -d $S/results ] && [ ! -d $S/results-610 ] && mv $S/results $S/results-610
mkdir -p $S/results $S/classII $S/mhc_fastas $S/immuannot_gtfs $S/figures $S/tables $S/workflow
rsync -a --exclude 'pggb_*/' $H/results-stage5/ $S/results/
rsync -a $H/results-stage5/pggb_* $S/results/ 2>/dev/null || true
rsync -a $H/classII5/classII.fa $H/classII5/classII.regions.tsv $H/classII5/out/ $S/classII/
rsync -a $H/mhc_all/ $S/mhc_fastas/
rsync -a $H/gtf_all/ $S/immuannot_gtfs/
rsync -a $H/viz5/figures/ $S/figures/
rsync -a $H/viz5/data/typing_concordance* $H/viz5/data/classII_haplotype_strings.tsv $H/viz5/data/classII_cluster_purity.tsv \
  $H/viz5/data/novel_alleles_table.tsv $H/viz5/data/novel_coding_all.tsv $H/viz5/data/hprc_r2_populations.tsv \
  $H/viz5/gene_copy_number.tsv $S/tables/
gzip -c $H/viz5/hla_calls.tsv > $S/tables/hla_calls.tsv.gz
rsync -a --exclude '.git' --exclude '__pycache__' $H/cwl/ $S/workflow/
cp $H/cwl/nig/inputs-stage5.yml $S/inputs-754-haplotypes.yml
chmod -R g+rX,o+rX $S 2>/dev/null || true
du -sh $S/* | sort -h | tail -12
echo SHAREDONE $(date)
