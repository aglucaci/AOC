#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------------------------
# AOC one-shot local runner (samples.csv-driven)
#
# What it does:
#   1) Ensures software/hyphy-analyses is present (clones if missing)
#   2) Runs the AOC workflow(s) using a samples.csv file (no yq config editing)
#
# Expected:
#   - Your Snakefile(s) read the samples sheet via config key: samples_csv
#   - You run from the repo root (or pass --workdir)
#
# Usage:
#   bash run_aoc.sh --samples samples.csv
#
# Optional:
#   --cores 8
#   --jobs 4
#   --snakefile workflow/Snakefile
#   --selection workflow/Snakefile_SelectionAnalysis
#   --summarize workflow/Snakefile_SummarizeResults
#   --latency-wait 60
#   --rerun-incomplete
#   --dry-run
# ------------------------------------------------------------------------------

# Pretty banner
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

# -----------------------------
# Defaults
# -----------------------------
SAMPLES_CSV=""
CORES="all"
JOBS="4"
LATENCY_WAIT="60"
KEEP_GOING="--keep-going"
REASON="--reason"
RERUN_INCOMPLETE="--rerun-incomplete"
DRYRUN=""

SNAKEFILE_MAIN="workflow/Snakefile"
SNAKEFILE_SELECTION="workflow/Snakefile_SelectionAnalysis"
SNAKEFILE_SUMMARIZE="workflow/Snakefile_SummarizeResults"

CLUSTER_CONFIG="config/cluster.json"   # kept for compatibility; safe to omit locally
USE_CLUSTER_CONFIG="yes"

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
    --selection)
      SNAKEFILE_SELECTION="$2"; shift 2;;
    --summarize)
      SNAKEFILE_SUMMARIZE="$2"; shift 2;;
    --no-cluster-config)
      USE_CLUSTER_CONFIG="no"; shift 1;;
    --workdir)
      WORKDIR="$2"; shift 2;;
    --no-rerun-incomplete)
      RERUN_INCOMPLETE=""; shift 1;;
    --dry-run|-n)
      DRYRUN="--dry-run"; shift 1;;
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
  --selection            Selection-analysis Snakefile (default: workflow/Snakefile_SelectionAnalysis)
  --summarize            Summarize-results Snakefile (default: workflow/Snakefile_SummarizeResults)
  --no-cluster-config    Do not pass --cluster-config config/cluster.json
  --workdir              Change into this directory before running
  --no-rerun-incomplete  Disable --rerun-incomplete
  --dry-run, -n          Snakemake dry run
  --help, -h             Show this help

Examples:
  bash scripts/run_aoc_samples.sh --samples samples.csv --cores 8 --jobs 4
  bash scripts/run_aoc_samples.sh -s config/samples.csv --no-cluster-config
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
# Ensure hyphy-analyses exists
# -----------------------------
#FOLDER="software/hyphy-analyses"
#URL="https://github.com/veg/hyphy-analyses.git"

#if [[ ! -d "${FOLDER}" ]]; then
#  echo "[INFO] Cloning hyphy-analyses into ${FOLDER} ..."
#  mkdir -p "$(dirname "${FOLDER}")"
#  git clone "${URL}" "${FOLDER}"
#else
#  echo "[INFO] Found ${FOLDER} (skipping clone)"
#fi

# -----------------------------
# Logs dir
# -----------------------------
# Handled internally in results/<sample>/logs
#mkdir -p logs

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

  echo "###############################################################################"
  echo "# ${title}"
  echo "###############################################################################"

  local cluster_args=()
  if [[ "${USE_CLUSTER_CONFIG}" == "yes" && -f "${CLUSTER_CONFIG}" ]]; then
    cluster_args=(--cluster-config "${CLUSTER_CONFIG}")
  fi

  # Always pass samples_csv into config for the Snakefile(s) to consume.
  snakemake \
    -s "${snakefile}" \
    "${cluster_args[@]}" \
    --jobs "${JOBS}" \
    --cores "${CORES}" \
    ${KEEP_GOING} \
    ${REASON} \
    --latency-wait "${LATENCY_WAIT}" \
    ${RERUN_INCOMPLETE} \
    ${DRYRUN} --printshellcmds \
    --config "samples_csv=${SAMPLES_CSV}"
  echo ""
}

# -----------------------------
# Run pipeline phases
# -----------------------------
run_smk "${SNAKEFILE_MAIN}"      "Running the AOC Snakemake pipeline (samples.csv)"
#run_smk "${SNAKEFILE_SELECTION}" "Selection analyses (recombination-free)"
#run_smk "${SNAKEFILE_SUMMARIZE}" "Visualization and summary"

echo "[DONE] All requested phases finished."
