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
# Analysis scripts (analysis/*.py): graph layouts, dendrogram, xlsx truth sets
OPENBLAS_NUM_THREADS=1 mm/envs/hla/bin/python3 -c "import networkx, Bio, scipy, openpyxl" 2>/dev/null || \
  OPENBLAS_NUM_THREADS=1 mm/envs/hla/bin/python3 -m pip install -q networkx biopython scipy openpyxl
# HLA typing env: FuFiHLA (long-read full-field typing, run on assembly pseudo-reads) and T1K (short reads)
[ -d mm/envs/typing ] || bin/micromamba create -y -r mm -n typing -c conda-forge -c bioconda fufihla t1k samtools
mkdir -p typing
[ -s typing/ref_data/ref.gene.fa.gz ] || (cd typing && ../mm/envs/typing/bin/fufihla-ref-prep)            # IPD-IMGT/HLA 3.65 at setup time
[ -s typing/t1k/hlaidx/hlaidx_dna_seq.fa ] || { mkdir -p typing/t1k && (cd typing/t1k && PATH=$H/mm/envs/typing/bin:$PATH t1k-build.pl -o hlaidx --download IPD-IMGT/HLA); }
# Minigraph-Cactus (whole-MHC pangenome graph for vg giraffe), static binary release
if [ ! -x cactus/cactus-bin-v3.3.0/bin/cactus_consolidated ]; then
  mkdir -p cactus && cd cactus
  curl -sfL -o cactus-bin-v3.3.0.tar.gz https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v3.3.0/cactus-bin-v3.3.0.tar.gz
  tar xzf cactus-bin-v3.3.0.tar.gz && cd cactus-bin-v3.3.0
  /usr/bin/python3 -m venv venv-cactus-v3.3.0 && . venv-cactus-v3.3.0/bin/activate
  python3 -m pip install -q -U setuptools pip wheel && python3 -m pip install -q -U -r ./toil-requirement.txt && python3 -m pip install -q -U .
  deactivate; cd $H
fi
echo "typing and graph tools ready"
