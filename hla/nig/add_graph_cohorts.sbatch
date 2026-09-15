#!/bin/bash
#SBATCH --partition=asianhla-c32 --account=asianhla-group
#SBATCH -c 16 --mem=150G -t 2-00:00:00 -J hla_stage5 -o add_graph_cohorts.%j.log
# Add the graph-only cohorts (K-PanRef Korean, CPC Chinese) and rerun everything with HPRC split by population:
#  1. MHC extraction + Immuannot on the MHC path segments from nig/kpanref_mhc.sbatch and nig/cpc_mhc.sbatch -> results-graphcohorts
#  2. copy their MHC fastas / GTFs into mhc_all / gtf_all
#  3. stage 5: aggregation, gene cuts, bundles, graphs on all haplotypes (HPRC cohorts from data/hprc_r2_populations.tsv)
#  4. class II span (DRA..DMA) + pgr-tk principal-bundle decomposition -> classII5/
set -euo pipefail
H=$HOME/hla; CWL=$H/cwl
export PATH=$H/toil/bin:$H/pgr-tk-bin/release:$H/mm/envs/hla/bin:$PATH OPENBLAS_NUM_THREADS=1
KF=$H/kpanref/fasta; CF=$H/cpc/fasta
# records were written as sample#hap#locus_offset by an earlier version of kpanref_mhc.sbatch; extract_mhc.py adds the prefix itself
for f in $KF/*.mhcctg.fa; do sed -i 's/^>[^#]*#[0-9]#/>/' $f; rm -f $f.fai; samtools faidx $f; done
COMMON="--reference /home/asianhla/data/HPRC_r2/fasta/GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz --immuannot-dir $H/Immuannot --immuannot-ref $H/Data-2024Feb02"
R1=$H/results-graphcohorts
if [ ! -s $R1/hla_calls.tsv ]; then
  python3 $CWL/scripts/make_inputs.py --dir KPanRef=$KF --dir CPC=$CF $COMMON --out $CWL/nig/inputs-graphcohorts.yml
  bash $CWL/nig/run_nig.sbatch $CWL/nig/inputs-graphcohorts.yml $R1
fi
cp $R1/*.mhc.fa $R1/*.mhc.tsv $H/mhc_all/
cp $R1/*.gtf.gz $H/gtf_all/
echo "stage 1-2 done $(date)"
python3 $CWL/scripts/make_inputs.py \
  --dir APR=/home/asianhla/data/upload/APR/assemblies --dir HPRC_r2=/home/asianhla/data/HPRC_r2/fasta \
  --dir JaSaPaGe=/home/asianhla/data/JaSaPaGe/assembly_clean/fasta --dir KPanRef=$KF --dir CPC=$CF \
  --populations $CWL/data/hprc_r2_populations.tsv $COMMON \
  --mhc-dir $H/mhc_all --gtf-dir $H/gtf_all --out $CWL/nig/inputs-stage5.yml 2>&1 | grep -v -e skipped -e duplicate
rm -rf $H/results-stage5; bash $CWL/nig/run_nig.sbatch $CWL/nig/inputs-stage5.yml $H/results-stage5
echo "stage 3 done $(date)"
C=$H/classII5; rm -rf $C $H/jobstore-classII5; mkdir -p $C/out; cd $C
python3 $CWL/scripts/extract_span.py --calls $H/results-stage5/hla_calls.tsv --fasta $H/results-stage5/MHC.fa \
  --from-gene HLA-DRA --to-gene HLA-DMA --name classII --flank 20000 --outdir $C
sed "s#/home/leechuck/hla/classII/classII.fa#$C/classII.fa#; s/threads: 4/threads: 16/" $H/classII/job.yml > job.yml
toil-cwl-runner --batchSystem single_machine --maxCores ${SLURM_CPUS_PER_TASK:-16} --jobStore $H/jobstore-classII5 --workDir ${SLURM_TMPDIR:-/tmp} \
  --clean onSuccess --no-container --preserve-environment PATH --bypass-file-store --outdir $C/out --logFile $C/toil.log \
  $CWL/tools/pgr_pbundle.cwl job.yml > $C/cwl-output.json
echo PIPEDONE $(date)
