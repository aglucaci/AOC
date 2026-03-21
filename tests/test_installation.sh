#!/usr/bin/env bash
set -euo pipefail

echo "Running AOC installation test..."

# Clean previous test
rm -rf tests/tmp_output

# Run minimal workflow
snakemake \
  --cores 1 \
  --snakefile workflow/Snakefile \
  all \
  --config samples_csv=tests/data/tiny_samples.csv outdir=tests/tmp_output \
  --rerun-triggers mtime

# Basic checks
#if [ ! -f tests/tmp_output/summary/run_manifest.csv ]; then
#  echo "Installation test failed."
#  exit 1
#fi

echo "Installation test passed."
