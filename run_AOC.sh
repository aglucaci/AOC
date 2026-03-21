#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# AOC one-shot local runner (samples.csv-driven)
#
# What it does:
#   1) Runs the AOC workflow(s) using a samples.csv file (no yq config editing)
#
# Expected:
#   - Your Snakefile(s) read the samples sheet via config key: samples_csv
#   - You run from the repo root (or pass --workdir)
#
# Usage:
#   bash run_AOC.sh --samples samples.csv
#
# Optional:
#   --cores 8
#   --jobs 4
#   --snakefile workflow/Snakefile
#   --latency-wait 60
#   --rerun-incomplete
#   --dry-run
#   --software-deployment-method conda
# ------------------------------------------------------------------------------

# Pretty banner
if [[ -t 1 && -z "${SLURM_JOB_ID:-}" ]]; then
clear || true
echo ""
cat <<'EOF'
     █████╗   ██████╗   ██████╗
    ██╔══██╗ ██╔═══██╗ ██╔════╝
    ███████║ ██║   ██║ ██║
    ██╔══██║ ██║   ██║ ██║
    ██║  ██║ ╚██████╔╝ ╚██████╔╝
    ╚═╝  ╚═╝  ╚═════╝   ╚═════╝
EOF
echo ""
fi

# -----------------------------
# Defaults
# -----------------------------
SAMPLES_CSV="samples.csv"
CORES="all"
JOBS="4"
LATENCY_WAIT="60"
KEEP_GOING="--keep-going"
REASON="--reason"
RERUN_INCOMPLETE="--rerun-incomplete"
DRYRUN=""
SOFTWARE_DEPLOYMENT_METHOD=""

SNAKEFILE_MAIN="workflow/Snakefile"

WORKDIR=""

export GRPC_VERBOSITY=NONE
export GRPC_TRACE=
export GLOG_minloglevel=3

# -----------------------------
# Args
# -----------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples|-s)
      SAMPLES_CSV="$2"; shift 2;;
    --cores|-c)
      CORES="$2"; shift 2;;
    --jobs|-j)
      JOBS="$2"; shift 2;;
    --latency-wait)
      LATENCY_WAIT="$2"; shift 2;;
    --snakefile)
      SNAKEFILE_MAIN="$2"; shift 2;;
    --workdir)
      WORKDIR="$2"; shift 2;;
    --no-rerun-incomplete)
      RERUN_INCOMPLETE=""; shift 1;;
    --dry-run|-n)
      DRYRUN="--dry-run"; shift 1;;
    --software-deployment-method)
      SOFTWARE_DEPLOYMENT_METHOD="--software-deployment-method $2"; shift 2;;
    --help|-h)
      cat <<EOF
Usage: bash scripts/run_aoc_samples.sh --samples samples.csv [options]

Required:
  --samples, -s          Path to samples.csv

Options:
  --cores, -c            Snakemake --cores value (default: all)
  --jobs, -j             Snakemake --jobs value (default: 4)
  --latency-wait         Snakemake --latency-wait (default: 60)
  --snakefile            Main Snakefile path (default: workflow/Snakefile)
  --workdir              Change into this directory before running
  --no-rerun-incomplete  Disable --rerun-incomplete
  --dry-run, -n          Snakemake dry run
  --software-deployment-method  Pass through Snakemake deployment backend(s), e.g. 'conda'
  --help, -h             Show this help

Examples:
  bash scripts/run_AOC.sh --samples samples.csv --cores 8 --jobs 4
  bash scripts/run_AOC.sh -s config/samples.csv
EOF
      exit 0;;
    *)
      echo "[ERROR] Unknown arg: $1" >&2
      echo "        Run with --help" >&2
      exit 2;;
  esac
done

if [[ -z "${SAMPLES_CSV}" ]]; then
  echo "[ERROR] --samples is required. Example: --samples samples.csv" >&2
  exit 2
fi

if [[ -n "${WORKDIR}" ]]; then
  cd "${WORKDIR}"
fi

if [[ ! -f "${SAMPLES_CSV}" ]]; then
  echo "[ERROR] samples.csv not found: ${SAMPLES_CSV}" >&2
  exit 2
fi

# -----------------------------
# Snakemake wrapper
# -----------------------------
run_smk () {
  local snakefile="$1"
  local title="$2"

  if [[ ! -f "${snakefile}" ]]; then
    echo "[WARN] Snakefile not found: ${snakefile} (skipping: ${title})"
    return 0
  fi

  #echo "###############################################################################"
  #echo "# ${title}"
  #echo "###############################################################################"

  # Always pass samples_csv into config for the Snakefile(s) to consume.
  snakemake \
    -s "${snakefile}" \
    --jobs "${JOBS}" \
    --cores "${CORES}" \
    ${KEEP_GOING} \
    --latency-wait "${LATENCY_WAIT}" \
    ${RERUN_INCOMPLETE} \
    ${SOFTWARE_DEPLOYMENT_METHOD} \
    all \
    ${DRYRUN} --printshellcmds \
    --config "samples_csv=${SAMPLES_CSV}" \
    --rerun-triggers mtime
  echo ""
}

# -----------------------------
# Run pipeline phases
# -----------------------------
run_smk "${SNAKEFILE_MAIN}"      "Running the AOC Snakemake workflow"

echo "[DONE] All requested phases finished."
