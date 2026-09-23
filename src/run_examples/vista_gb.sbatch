#!/bin/bash
#SBATCH -J gschism          # job name
#SBATCH -p gb               # GH200 queue
#SBATCH -N 1                # nodes
#SBATCH -n 1                # tasks
##SBATCH --gres=gpu:2        # gb already provides the GPU; uncomment if nvidia-smi is empty
#SBATCH -t 00:20:00         # wall time
#SBATCH -o %x.%j.out        # stdout
#SBATCH -e %x.%j.err        # stderr
#
# Usage:
#   sbatch vista_gb.sbatch CASEDIR
#   sbatch vista_gb.sbatch CASEDIR OUTDIR
#
# CASEDIR is the SCHISM input directory (hgrid.gr3, param.nml, ...).
# OUTDIR defaults to CASEDIR/gpuoutput if omitted.

set -euo pipefail

# Absolute path on VISTA (do not derive from $0; sbatch runs a spool copy).
SRC_DIR="/home1/10922/jmglez/src/gschism/src"
ENV_NAME="cupy311"
MAMBA_EXE="$HOME/bin/micromamba"
ENV_PREFIX="$HOME/.local/share/mamba/envs/${ENV_NAME}"

if [[ $# -lt 1 ]]; then
    echo "usage: sbatch $0 CASEDIR [OUTDIR]" >&2
    exit 1
fi

CASE_DIR="$(readlink -f "$1")"
if [[ $# -ge 2 ]]; then
    OUT_DIR="$(readlink -f "$2")"
else
    OUT_DIR=""
fi

if [[ ! -x "$MAMBA_EXE" ]]; then
    echo "micromamba not found: $MAMBA_EXE" >&2
    exit 1
fi

if [[ ! -d "$ENV_PREFIX" ]]; then
    echo "env '${ENV_NAME}' not found: $ENV_PREFIX" >&2
    echo "create it with: micromamba env create -n ${ENV_NAME} -f ${SRC_DIR}/environment.yaml" >&2
    exit 1
fi

if [[ ! -d "$CASE_DIR" ]]; then
    echo "case directory not found: $CASE_DIR" >&2
    exit 1
fi

module purge
module load TACC
module load gcc/15.1.0
module load cuda/13.1
module load metis || true

export CC=$(which gcc)
export CXX=$(which g++)

echo "TACC_CUDA_DIR=$TACC_CUDA_DIR"
nvcc --version || true
nvidia-smi || true

export MAMBA_EXE
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba activate "$ENV_NAME"

# Match math_libs to the loaded CUDA; stub dirs break cuBLAS at runtime
export CUDA_HOME="$TACC_CUDA_DIR"
CUDA_VER="$(basename "$TACC_CUDA_DIR")"
NVIDIA_ROOT="$(dirname "$(dirname "$TACC_CUDA_DIR")")"
export MATH_LIB="$NVIDIA_ROOT/math_libs/${CUDA_VER}/lib64"
export LD_LIBRARY_PATH="$MATH_LIB:$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$(echo "$LD_LIBRARY_PATH" | tr ':' '\n' | grep -v '/lib64/stubs' | paste -sd:)"

export PYTHONPATH="${SRC_DIR}:${PYTHONPATH:-}"
cd "$SRC_DIR"

echo "Starting at: $(date)"
echo "Host: $(hostname)"
echo "CASE_DIR=$CASE_DIR"
echo "OUT_DIR=${OUT_DIR:-$CASE_DIR/gpuoutput}"

PY="$SRC_DIR/python/main.py"
if [[ -n "$OUT_DIR" ]]; then
    srun python -u "$PY" "$CASE_DIR" --outDir "$OUT_DIR"
else
    srun python -u "$PY" "$CASE_DIR"
fi

echo "Done at: $(date)"
