
<p align="center">
  <img src="logo/AOC_Logo_3.png" alt="AOC Logo" width="400"/>
</p>

# Analysis of Orthologous Collections (AOC)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://github.com/aglucaci/Analysis-of-Orthologous-Collections/blob/main/LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/aglucaci/Analysis-of-Orthologous-Collections)](https://github.com/aglucaci/Analysis-of-Orthologous-Collections/commits/main)
[![GitHub issues](https://img.shields.io/github/issues/aglucaci/Analysis-of-Orthologous-Collections)](https://github.com/aglucaci/Analysis-of-Orthologous-Collections/issues)
[![Snakemake](https://img.shields.io/badge/Snakemake-pipeline-brightgreen)](https://snakemake.readthedocs.io/)
[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)


**AOC** is a reproducible, Snakemake-based pipeline for the automated investigation of molecular evolution across orthologous protein-coding genes. It integrates high-throughput alignment, recombination detection, phylogenetic reconstruction, and a comprehensive suite of HyPhy-based selection analyses.

---

## Features

- Codon-aware multiple sequence alignment using **MACSE2**
- Recombination detection with **HyPhy GARD**
- Phylogenetic tree inference using **IQ-TREE**
- Comprehensive molecular evolution analysis via **HyPhy**:
  - Site-level: MEME, FEL, SLAC, FUBAR
  - Branch-level: aBSREL, BUSTED
  - Gene-wide: Model testing, RELAX
  - Co-evolution and heterogeneity: BGM, FMM
- Lineage assignment and tree annotation via **NCBI Taxonomy + ete3**
- Automated result summarization and visualization
- Compatible with local or HPC environments

---

## Directory Structure

```
AOC/
├── config/                # Configuration YAMLs and environment files
├── data/                  # Example input datasets (e.g., Primate_ACE2/)
├── results/               # Output directory for all results
├── scripts/               # Custom helper scripts (taxonomy, annotation, summaries)
├── workflow/              # Snakemake rules and pipeline logic
├── Launch_AOC_Locally.sh  # Run script for local machines
├── Launch_AOC_HPC.sh      # Run script for SLURM clusters
└── README.md
```

---

## Installation

We recommend using **conda** for environment management.

```bash
conda env create -f config/environment.yml
conda activate AOC
```

---

## Input Requirements

Each dataset should include:
- A protein FASTA file of orthologs
- A matching transcript FASTA file
- A metadata CSV file with RefSeq accessions

These files are downloaded from the NCBI Orthologs Database, which provides curated orthologous gene information across species. For example, orthologs of the human ACE2 gene (Gene ID: 59272) within primates (TaxID: 9443) can be retrieved via this interface.

Full search link for **Primate ACE2** example: [https://www.ncbi.nlm.nih.gov/gene/59272/ortholog/?scope=9443&term=ACE2](https://www.ncbi.nlm.nih.gov/gene/59272/ortholog/?scope=9443&term=ACE2)


Example:

```
data/Primate_ACE2/
├── ACE2_orthologs.csv
├── ACE2_refseq_protein.fasta
└── ACE2_refseq_transcript.fasta
```

---


## Configuring Your Analysis: Editing the YAML File
Before running the workflow, you must specify the genes you want to analyze by editing the YAML configuration file (`user/genes.yml`).

Navigate to the section labeled:

```
# =============================================================================
# Multiple Genes
# =============================================================================

# Edit the 'GENES' variable below,
#     these should be the names of folders in the 'data' directory

# -----------------------------------------------------------------------------
# For multiple GENES: use a comma-delimited list
# GENES: Primate_ACE2,Primate_TP53,Primate_BTG1,Primate_REM2

# If you only want to run one gene.
#
GENES: Primate_ACE2
```

### Instructions:
Each entry in the `GENES:` line should match the name of a folder inside the (`data/ directory`).

These folders must contain:

* A protein FASTA file
* A transcript FASTA file
* A metadata CSV file

To analyze multiple genes, provide a comma-separated list:

`GENES: Primate_ACE2,Primate_TP53,Primate_BTG1`

To analyze only one gene:

`GENES: Primate_ACE2`

**Reminder!** Each gene folder should correspond to ortholog datasets downloaded from the [NCBI Orthologs Database](https://www.ncbi.nlm.nih.gov/gene/59272/ortholog/?scope=7742&term=ACE2).

## Running the Pipeline

Before running any scripts, ensure you are in the root directory of the repository.

### Local Execution (Multiple genes at a time)

```
bash Launch_AOC_Locally.sh
```

### HPC Execution (SLURM)

The configuration details described are specific for a SLURM-based HPC system. This general framework will need to be modified for other types of schedulers (PBS/Torque, SGE, HTCondor, Kubernetes).

#### HPC Setup – Required (`config/cluster.json` file
If you're planning to run the pipeline on a high-performance computing (HPC) cluster, you must provide an updated ('config/cluster.json`) file. This file defines the default resource allocation for each job submitted by Snakemake.

A minimal working example looks like this:

```

{
  "__default__": {
    "cluster": "sbatch",
    "nodes": 1,
    "ppn": 8,
    "name": "scu-cpu",
    "walltime": "72:00:00"
  }
}
```

**What Each Field Means:**

* (`"cluster": "sbatch"`) — Tells Snakemake to use SLURM (sbatch) for job submission.
* (`"nodes": 1`) — Number of nodes to allocate.
* (`"ppn": 8`) — Processors per node (can also be cpus-per-task depending on SLURM setup).
* (`"name": "scu-cpu"`) — Job name prefix; customize it for easier tracking.
* (`"walltime": "72:00:00"`) — Maximum run time (in HH:MM:SS) for each job.


After this is configured Launch the Snakemake via: 
```
sbatch Launch_AOC_HPC.sh
```

## Checklist Before Running
*  You’ve edited (`user/genes.yml`) with the correct `GENES:` variable
*  Each gene listed exists as a folder in the (`data/`) directory
*  Each folder contains protein FASTA, transcript FASTA, and metadata CSV
*  You're in the correct working directory
*  Your environment is activated (conda activate AOC, etc.)

---

## Output

All results are saved to `results/<GENE>/` and include:
- Codon alignments
- Phylogenetic trees (`.treefile`)
- JSON and CSV summaries from each HyPhy method
- Annotated trees with NCBI taxonomic metadata
- Visual summaries of sites under selection (`Visualization/`)
- Summary statistics and results in CSV-Format (`Tables/`)

---

## Methods Summary

| Method   | Purpose                                  | Scale         |
|----------|------------------------------------------|----------------|
| FEL      | Pervasive site-level selection (ML)      | Site           |
| SLAC     | Fast site-level selection (counting)     | Site           |
| MEME     | Episodic (branch-specific) selection     | Site           |
| FUBAR    | Bayesian site-level selection            | Site           |
| aBSREL   | Adaptive branch-level selection          | Branch         |
| BUSTED   | Gene-wide episodic positive selection    | Gene/Branch    |
| RELAX    | Tests relaxation/intensification of ω    | Lineage        |
| BGM      | Detects co-evolving sites                | Site-pair      |
| FMM      | Finite mixture model of site classes     | Site           |

---

## Example Use Case

We include a case study on **Primate ACE2** evolution:

This generates:

* Site-level selection maps
* Annotated phylogenetic trees
* Tables of positively/negatively selected sites across species
* Lineage-specific selection comparisons

---

## Citation

If you use AOC in your work, please cite:

> Lucaci AG, Pond SLK. AOC: Analysis of Orthologous Collections - an application for the characterization of natural selection in protein-coding sequences. ArXiv [Preprint]. 2024 Jun 13:arXiv:2406.09522v1. PMID: 38947939; PMCID: PMC11213150.

---

## Contact

Created and maintained by **Alexander G. Lucaci**  
Questions? Feature requests? Open an [issue](https://github.com/aglucaci/Analysis-of-Orthologous-Collections/issues) or contact [agl4001@med.cornell.edu](mailto:agl4001@med.cornell.edu)

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

---
