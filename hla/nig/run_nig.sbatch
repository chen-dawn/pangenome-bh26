#!/bin/bash
#SBATCH --job-name=hla_pangenome
#SBATCH --partition=asianhla-c32
#SBATCH --account=asianhla-group
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
# Run the HLA CWL workflow on the NIG BioHackathon node with Toil.
# The compute node has no container runtime, so tools come from the micromamba
# env built with nig/setup_env.sh and the runner uses --no-container.
# cwltool --parallel deadlocks on the 610-way scatter (stalled after ~390 jobs),
# hence toil-cwl-runner (single-machine batch system, restartable job store).
# usage: sbatch run_nig.sbatch <job.yml> <outdir> [restart]
set -euo pipefail
WF=${WF:-$HOME/hla/cwl}
JOB=${1:-$WF/nig/inputs-nig.yml}
OUT=${2:-$HOME/hla/results}
MODE=${3:-}
export PATH=$HOME/hla/toil/bin:$HOME/hla/pgr-tk-bin/release:$HOME/hla/mm/envs/hla/bin:$PATH
export TMPDIR=${SLURM_TMPDIR:-/tmp}
CORES=${SLURM_CPUS_PER_TASK:-28}
MEM=${SLURM_MEM_PER_NODE:-200000}   # MB
JOBSTORE=$HOME/hla/jobstore-$(basename "$OUT")
mkdir -p "$OUT"
echo "Job ${SLURM_JOB_ID:-local} on $(hostname), $(date), cores=$CORES mem=${MEM}M"
RESTART=""
if [ "$MODE" = "restart" ]; then RESTART="--restart"; else rm -rf "$JOBSTORE"; fi
toil-cwl-runner $RESTART \
  --batchSystem single_machine --maxCores "$CORES" --maxMemory "${MEM}M" --maxDisk 150G \
  --jobStore "$JOBSTORE" --workDir "$TMPDIR" --clean onSuccess \
  --no-container --preserve-environment PATH --retryCount 1 \
  --bypass-file-store --outdir "$OUT" --logLevel INFO --logFile "$OUT/toil.log" \
  "$WF/hla_pangenome.cwl" "$JOB" > "$OUT/cwl-output.json"
echo "done $(date)"
