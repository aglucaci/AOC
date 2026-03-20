#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# AOC stable installer
# - Avoids broken "mamba" installs (common on older Anaconda bases)
# - Supports explicit override via FRONTEND_OVERRIDE=conda|micromamba|mamba
# - Falls back safely to python-only dependencies if environment solve fails
# - Bootstraps hyphy-analyses repo if missing
#
# Usage:
#   bash ./install.sh AOC envs/AOC.yaml
#   FRONTEND_OVERRIDE=conda bash ./install.sh AOC envs/AOC.yaml
#   FRONTEND_OVERRIDE=micromamba bash ./install.sh AOC envs/AOC.yaml
# -----------------------------------------------------------------------------

ENV_NAME="${1:-AOC}"
ENV_FILE="${2:-envs/AOC.yaml}"

echo "[install] env_name=${ENV_NAME}"
echo "[install] env_file=${ENV_FILE}"

have() { command -v "$1" >/dev/null 2>&1; }

# --- sanity checks for frontends (important: mamba can exist but be broken)
mamba_ok() {
  have mamba || return 1
  mamba --version >/dev/null 2>&1 || return 1
  return 0
}
conda_ok() {
  have conda || return 1
  conda --version >/dev/null 2>&1 || return 1
  return 0
}
micromamba_ok() {
  have micromamba || return 1
  micromamba --version >/dev/null 2>&1 || return 1
  return 0
}

# --- choose frontend safely
FRONTEND_OVERRIDE="${FRONTEND_OVERRIDE:-}"
FRONTEND=""

if [ -n "${FRONTEND_OVERRIDE}" ]; then
  FRONTEND="${FRONTEND_OVERRIDE}"
  echo "[install] FRONTEND_OVERRIDE=${FRONTEND}"
else
  # Prefer micromamba (self-contained), then conda, then mamba (only if sane)
  if micromamba_ok; then
    FRONTEND="micromamba"
  elif conda_ok; then
    FRONTEND="conda"
  elif mamba_ok; then
    FRONTEND="mamba"
  else
    echo "[install] ERROR: No working conda frontend found (micromamba/conda/mamba)."
    echo "  Install one of:"
    echo "    - micromamba (recommended): https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
    echo "    - Miniforge/Miniconda"
    exit 1
  fi
fi

# Validate override
case "${FRONTEND}" in
  micromamba) micromamba_ok || { echo "[install] ERROR: micromamba not usable."; exit 1; } ;;
  conda)      conda_ok || { echo "[install] ERROR: conda not usable."; exit 1; } ;;
  mamba)      mamba_ok || { echo "[install] ERROR: mamba not usable (likely broken). Try FRONTEND_OVERRIDE=conda"; exit 1; } ;;
  *)          echo "[install] ERROR: unknown FRONTEND=${FRONTEND} (use conda|micromamba|mamba)"; exit 1 ;;
esac

echo "[install] using frontend: ${FRONTEND}"

# --- helpers to run commands in env
run_in_env() {
  # $1.. = command
  if [ "${FRONTEND}" = "micromamba" ]; then
    micromamba run -n "${ENV_NAME}" "$@"
  else
    conda run -n "${ENV_NAME}" "$@"
  fi
}

# --- create env (best effort)
create_env() {
  if [ "${FRONTEND}" = "micromamba" ]; then
    micromamba create -y -n "${ENV_NAME}" -f "${ENV_FILE}"
  elif [ "${FRONTEND}" = "mamba" ]; then
    # mamba uses conda-style syntax
    mamba env create -n "${ENV_NAME}" -f "${ENV_FILE}"
  else
    conda env create -n "${ENV_NAME}" -f "${ENV_FILE}"
  fi
}

update_env() {
  if [ "${FRONTEND}" = "micromamba" ]; then
    micromamba update -y -n "${ENV_NAME}" -f "${ENV_FILE}"
  elif [ "${FRONTEND}" = "mamba" ]; then
    mamba env update -n "${ENV_NAME}" -f "${ENV_FILE}"
  else
    conda env update -n "${ENV_NAME}" -f "${ENV_FILE}"
  fi
}

create_python_only_fallback() {
  echo "[install] Fallback: creating python-only environment (no HyPhy/IQ-TREE/FastTree)."
  if [ "${FRONTEND}" = "micromamba" ]; then
    micromamba create -y -n "${ENV_NAME}" python=3.11 pip
    micromamba run -n "${ENV_NAME}" pip install -U \
      snakemake pandas numpy biopython matplotlib pytest pyyaml rich
  else
    conda create -y -n "${ENV_NAME}" python=3.11 pip
    conda run -n "${ENV_NAME}" pip install -U \
      snakemake pandas numpy biopython matplotlib pytest pyyaml rich
  fi
  echo "[install] WARNING: Tool binaries (hyphy/iqtree/fasttree) not installed in fallback."
  echo "[install] If possible, install via conda-forge/bioconda later."
}

echo "[install] creating/updating environment…"

set +e
create_env
RC=$?
set -e

if [ $RC -ne 0 ]; then
  echo "[install] env create failed (rc=${RC}); trying env update…"
  set +e
  update_env
  RC2=$?
  set -e
  if [ $RC2 -ne 0 ]; then
    echo "[install] env update also failed (rc=${RC2}); using python-only fallback…"
    create_python_only_fallback
  fi
fi

echo "[install] environment ready: ${ENV_NAME}"

# --- install hyphy-analyses (bf files)
#if [ ! -d "software/hyphy-analyses" ]; then
#  echo "[install] hyphy-analyses not found; cloning…"
#  mkdir -p software
#  if have git; then
#    set +e
#    git clone --depth 1 https://github.com/veg/hyphy-analyses.git software/hyphy-analyses
#    GRC=$?
#    set -e
#    if [ $GRC -ne 0 ]; then
#      echo "[install] WARNING: failed to clone hyphy-analyses."
#      echo "  You can manually clone to: software/hyphy-analyses"
#    fi
#  else
#    echo "[install] WARNING: git not found; cannot clone hyphy-analyses automatically."
#    echo "  Please install git or manually download: software/hyphy-analyses"
#  fi
#else
#  echo "[install] hyphy-analyses present."
#fi

# --- quick smoke checks (non-fatal)
echo "[install] smoke checks…"
set +e
run_in_env python -V
run_in_env snakemake --version
run_in_env pytest --version
set -e

echo "[install] done."
