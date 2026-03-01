<p align="center">
  <img src="https://raw.githubusercontent.com/aglucaci/AOC/refs/heads/develop/logo/aoc_logo.png" alt="AOC" width="400"/>
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
├── workflow/
│   ├── Snakefile
│   ├── Snakefile_SelectionAnalysis
│   └── Snakefile_SummarizeResults
├── config/
│   └── config.yaml
├── scripts/
├── tests/
├── envs/
│   └── AOC.yaml
├── install.sh
├── run_aoc.sh
├── submit_aoc.slurm
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

**On MAC OSX (do this first)**
```
conda config --env --set subdir osx-64
```

**Run installation script**

```
bash install.sh aoc envs/aoc.yaml
```

### Optional: Force a frontend

```bash
FRONTEND_OVERRIDE=conda bash install.sh aoc envs/aoc.yaml
FRONTEND_OVERRIDE=micromamba bash install.sh aoc envs/aoc.yaml
```

After installation:

```
conda activate aoc
```

---

## Input Format

AOC is driven by a `samples.csv` file.

### Required columns

```
sample,codon_fasta
```

### Example

```
BDNF,data/BDNF.fasta
TP53,data/TP53.fasta
```

Each row corresponds to one ortholog dataset.

Optionally add a column: `sequence_labels_csv`

```
sample,codon_fasta,sequence_labels_csv
BDNF-annot,data/BDNF/BDNF.small.fasta,data/BDNF/BDNF.small.sequence_labels.csv
```

And format the `sequence_labels_csv` file in this format, with `label` corresponding the branch label, and `fasta_sequence_header` corresponding to the fasta sequence header description:

```
label,fasta_sequence_header
Test,"NM_001709.5 Homo sapiens brain derived neurotrophic factor (BDNF), transcript variant 4, mRNA"
Test,"NM_001270630.1 Rattus norvegicus brain-derived neurotrophic factor (Bdnf), transcript variant 1, mRNA"
Background,"XM_011226480.3 PREDICTED: Ailuropoda melanoleuca brain derived neurotrophic factor (BDNF), mRNA"
Background,"XM_007497196.2 PREDICTED: Monodelphis domestica brain-derived neurotrophic factor (BDNF), transcript variant X1, mRNA"
Background,"NM_001081787.1 Equus caballus brain derived neurotrophic factor (BDNF), mRNA"
```

---

## Running the Pipeline

From the root directory:

```bash
bash run_aoc.sh --samples samples.csv
```

---

## HPC Execution (SLURM Example)

```bash
sbatch submit_aoc.slurm
```

---

## Selection Methods Included

| Method       | Scale   | Purpose |
|--------------|----------|----------|
| FEL          | Site     | Pervasive selection |
| MEME         | Site     | Episodic selection |
| aBSREL       | Branch   | Adaptive branch selection |
| BUSTED-S-MH  | Gene     | Gene-wide episodic selection |
| CFEL         | Site     | Contrast site-level selection |
| RELAX        | Lineage  | Selection intensity shifts |

---

## Automated Testing

Quick installation validation:

```
bash tests/test_installation.sh
```

---

## Citation

If you use AOC in your work, please cite:

```
Lucaci AG, Pond SLK. AOC: Analysis of Orthologous Collections - an application for the characterization of natural selection in protein-coding sequences. ArXiv [Preprint]. 2024 Jun 13:arXiv:2406.09522v1. PMID: 38947939; PMCID: PMC11213150.
```
---

## License

GPL-3.0
