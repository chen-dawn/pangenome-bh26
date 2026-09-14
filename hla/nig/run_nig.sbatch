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
# Run the HLA CWL workflow on the NIG BioHackathon node.
# The compute node has no container runtime, so tools come from the micromamba
# env built with nig/setup_env.sh and cwltool runs with --no-container.
set -euo pipefail
WF=${WF:-$HOME/hla/cwl}
JOB=${1:-$WF/nig/inputs-nig.yml}
OUT=${2:-$HOME/hla/results}
export PATH=$HOME/hla/pgr-tk-bin/release:$HOME/hla/mm/envs/hla/bin:$PATH
export TMPDIR=${TMPDIR:-$HOME/hla/tmp}
mkdir -p "$OUT" "$TMPDIR"
echo "Job ${SLURM_JOB_ID:-local} on $(hostname), $(date)"
cwltool --parallel --no-container --preserve-environment PATH \
  --tmpdir-prefix "$TMPDIR/" --tmp-outdir-prefix "$TMPDIR/" \
  --outdir "$OUT" "$WF/hla_pangenome.cwl" "$JOB" > "$OUT/cwl-output.json" 2> "$OUT/cwltool.log"
echo "done $(date)"
