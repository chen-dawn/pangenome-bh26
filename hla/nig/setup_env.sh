#!/bin/bash
# One-time user-space environment on the NIG login node (no container runtime on the compute node).
set -euo pipefail
H=$HOME/hla; mkdir -p "$H"; cd "$H"
[ -x bin/micromamba ] || curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
export MAMBA_ROOT_PREFIX=$H/mm
[ -d mm/envs/hla ] || bin/micromamba create -y -r mm -n hla -c conda-forge -c bioconda \
  minimap2=2.28 samtools=1.21 seqkit pggb=0.7.4 odgi vg wfmash bedtools gfatools \
  python=3.11 pandas matplotlib biopython pyyaml
[ -d Immuannot ] || git clone https://github.com/YingZhou001/Immuannot.git
[ -d Data-2024Feb02 ] || { curl -sL -o Data-2024Feb02.tar.gz https://zenodo.org/api/records/10948964/files/Data-2024Feb02.tar.gz/content; tar xzf Data-2024Feb02.tar.gz; }
echo "env ready: $H"
# pgr-tk 0.5.1 command-line binaries (bioconda ships only the Python module)
[ -x pgr-tk-bin/release/pgr-query ] || { curl -sL -o pgr-tk-v0.5.1.zip https://github.com/GeneDx/pgr-tk/releases/download/v0.5.1/pgr-tk-v0.5.1.zip; mkdir -p pgr-tk-bin; (cd pgr-tk-bin && unzip -qo ../pgr-tk-v0.5.1.zip && chmod +x release/*); }
# Toil CWL runner (cwltool --parallel deadlocks on large scatters)
[ -x toil/bin/toil-cwl-runner ] || { python3 -m venv toil && toil/bin/pip install -q --upgrade pip && toil/bin/pip install -q "toil[cwl]"; }
# The Toil venv is first on PATH in run_nig.sbatch, so its python3 runs the
# aggregation script (--no-container): it needs pandas + matplotlib too.
toil/bin/python3 -c "import pandas, matplotlib" 2>/dev/null || toil/bin/pip install -q pandas matplotlib
