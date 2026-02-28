<p align="center">
  <img src="logo/AOC_Logo_3.png" alt="AOC Logo" width="400"/>
</p>

# Analysis of Orthologous Collections (AOC)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)]
[![Snakemake](https://img.shields.io/badge/Snakemake-pipeline-brightgreen)]
[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)]

---

## Overview

**AOC** is a reproducible, modular Snakemake workflow for ortholog-aware evolutionary analysis of protein-coding genes.

AOC integrates:

- Codon-aware alignment
- Phylogenetic reconstruction
- Recombination detection
- Branch labeling
- HyPhy-based molecular evolution analyses
- Automated summarization and reporting

The workflow is designed for both local execution and HPC environments, and scales across many ortholog datasets via a `samples.csv` configuration.

---

## Core Features

- Codon-aware alignments (MACSE2)
- Phylogenetic inference (IQ-TREE)
- Recombination detection (HyPhy GARD)
- Comprehensive selection inference:
  - MG94
  - FEL
  - MEME
  - CFEL
  - RELAX
  - aBSREL
  - BUSTED-S-MH
- Branch labeling workflows
- JSON parsing and summarization
- Automated run manifest generation
- Pytest-based automated testing
- Stable environment installer (`install.sh`)

---

## Repository Structure

```
AOC/
├── Snakefile
├── workflow/
│   ├── Snakefile_SelectionAnalysis
│   └── Snakefile_SummarizeResults
├── config/
│   └── config.yaml
├── scripts/
├── tests/
├── envs/
│   └── AOC.yaml
├── install.sh
└── README.md
```

---

## Installation (Recommended)

AOC provides a stable installer that:

- Detects micromamba / conda / mamba safely
- Avoids broken mamba installations
- Creates or updates environments
- Falls back to Python-only mode if solver fails
- Performs smoke testing

### Run installer

```bash
bash install.sh aoc envs/AOC.yaml
```

### Optional: Force a frontend

```bash
FRONTEND_OVERRIDE=conda bash install.sh aoc envs/AOC.yaml
FRONTEND_OVERRIDE=micromamba bash install.sh aoc envs/AOC.yaml
```

After installation:

```bash
conda activate aoc
```

---

## Input Format

AOC is driven by a `samples.csv` file.

### Required columns

```
sample,fasta
```

### Example

```
BDNF,data/BDNF.fasta
TP53,data/TP53.fasta
```

Each row corresponds to one ortholog dataset.

---

## Running the Pipeline

From the root directory:

```bash
snakemake --cores 8   --configfile config/config.yaml   --config samples_csv=samples.csv   --rerun-incomplete   --printshellcmds
```

### Override output directory

```bash
--config outdir=results
```

---

## HPC Execution (SLURM Example)

```bash
snakemake   --jobs 100   --cluster "sbatch --cpus-per-task={threads} --time=48:00:00"   --configfile config/config.yaml   --config samples_csv=samples.csv
```

---

## Output Structure

```
<OUTDIR>/
├── alignments/
├── trees/
├── selection/<sample>/
├── labels/
├── summary/
└── logs/
```

Key outputs include:

- Codon alignments
- Phylogenetic trees (.nwk)
- Selection model JSON outputs
- Summarized CSV tables
- run_manifest.csv
- Log files

---

## Selection Methods Included

| Method       | Scale   | Purpose |
|--------------|----------|----------|
| MG94         | Gene     | Baseline codon model |
| FEL          | Site     | Pervasive selection |
| MEME         | Site     | Episodic selection |
| CFEL         | Site     | Contrast site-level selection |
| RELAX        | Lineage  | Selection intensity shifts |
| aBSREL       | Branch   | Adaptive branch selection |
| BUSTED-S-MH  | Gene     | Gene-wide episodic selection |

---

## Automated Testing

Run full test suite:

```bash
pytest
```

Quick installation validation:

```bash
bash scripts/test_installation.sh
```

---

## Citation

If you use AOC in your work, please cite:

Lucaci AG, Pond SLK. AOC: Analysis of Orthologous Collections. 2024.

---

## License

GPL-3.0
