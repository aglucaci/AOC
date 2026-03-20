<p align="center">
  <img src="https://raw.githubusercontent.com/aglucaci/AOC/refs/heads/develop/logo/AOC_logo.png" alt="AOC" width="400"/>
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

## Statement of Need
While individual HyPhy analyses can be run through DataMonkey or the HyPhy command line, AOC is designed for reproducible large-scale analyses across many genes or datasets. It automates alignment preparation, phylogenetic inference, branch labeling, multiple selection tests, and standardized result aggregation within a single workflow.

---

## Core Features

- Codon-aware alignments (MACSE2)
- Phylogenetic inference (IQ-TREE)
- Recombination detection (HyPhy GARD)
- Comprehensive selection inference:
  - FEL
  - MEME
  - CFEL
  - RELAX
  - aBSREL
  - BUSTED-S-MH
- Branch labeling workflows
- JSON parsing and summarization
- Automated run manifest generation
- Stable environment installer (`install.sh`)

---

## Repository Structure

```
AOC/
├── workflow/
│   ├── Snakefile
├── config/
│   └── config.yaml
├── scripts/
├── tests/
├── envs/
│   └── AOC.yaml
├── install.sh
├── run_AOC.sh
├── submit_AOC.slurm
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
bash install.sh AOC envs/AOC.yaml
```

### Optional: Force a frontend

```bash
FRONTEND_OVERRIDE=conda bash install.sh AOC envs/AOC.yaml
FRONTEND_OVERRIDE=micromamba bash install.sh AOC envs/AOC.yaml
```

After installation:

```
conda activate AOC
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

Branches labeled “Test” represent the foreground lineages where a specific evolutionary hypothesis (e.g., adaptive selection) is being evaluated, while “Background” branches represent the remainder of the phylogeny and serve as a reference group against which evolutionary patterns in the Test set are compared.

---

## Running the Pipeline

From the root directory:

```bash
bash run_AOC.sh --samples samples.csv
```

---

## HPC Execution (SLURM Example)

```bash
sbatch submit_AOC.slurm
```

The example `submit_AOC.slurm` writes Slurm stdout/stderr to `AOC_<jobid>.out` and `AOC_<jobid>.err` in the submission directory, so no pre-existing `logs/` folder is required.

---

## Selection Methods Included

AOC integrates several widely used codon-based evolutionary models implemented in
[HyPhy (Hypothesis Testing using Phylogenies)](https://hyphy.org) to detect signals
of natural selection across protein-coding genes. These approaches operate at
different biological scales (site, branch, lineage, and gene) and capture
complementary evolutionary signals. Using multiple tests together improves
robustness because positive selection can manifest differently depending on the
evolutionary scenario.

All models are based on codon substitution frameworks (typically variants of
MG94 or related codon models) that estimate the ratio of nonsynonymous to
synonymous substitution rates (**dN/dS**, also called **ω**).

Interpretation of ω:

- **ω > 1** → positive (diversifying) selection  
- **ω = 1** → neutral evolution  
- **ω < 1** → purifying selection  

Each HyPhy method tests a different hypothesis about how selection acts across
sites and lineages.

| Method | Scale | Purpose |
|------|------|------|
| FEL | Site | Pervasive selection |
| MEME | Site | Episodic selection |
| aBSREL | Branch | Adaptive branch selection |
| BUSTED-S-MH | Gene | Gene-wide episodic selection |
| CFEL | Site | Contrast site-level selection |
| RELAX | Lineage | Selection intensity shifts |

---

### FEL — Fixed Effects Likelihood

**FEL (Fixed Effects Likelihood)** tests for **pervasive selection at individual
codon sites**. It estimates synonymous and nonsynonymous substitution rates
independently for each site using maximum likelihood.

Key characteristics:

- Detects **consistent selection across the entire phylogeny**
- Identifies sites under **persistent positive or negative selection**
- Conservative but interpretable site-level estimates

FEL is most appropriate when the selective pressure is expected to be **stable
across evolutionary time**.

Documentation:  
https://hyphy.org/methods/selection-methods/#fel

---

### MEME — Mixed Effects Model of Evolution

**MEME (Mixed Effects Model of Evolution)** detects **episodic positive selection
at individual sites**. Unlike FEL, MEME allows selection to occur **only on a
subset of branches**.

Key characteristics:

- Detects **transient or lineage-specific adaptive events**
- Combines site-level and branch-level modeling
- Powerful for detecting **adaptive bursts**

MEME is widely used when adaptive events are expected to occur **sporadically
during evolution**.

Documentation:  
https://hyphy.org/methods/selection-methods/#meme

---

### aBSREL — Adaptive Branch-Site Random Effects Likelihood

**aBSREL** identifies **branches experiencing episodic diversifying selection**
across a gene.

Key characteristics:

- Tests each branch independently
- Allows multiple ω rate classes on each branch
- Detects **adaptive episodes affecting subsets of sites**

This method is useful for identifying **specific evolutionary lineages
undergoing adaptation**.

Documentation:  
https://hyphy.org/methods/selection-methods/#absrel

---

### BUSTED-S-MH — Branch-Site Unrestricted Statistical Test for Episodic Diversification

**BUSTED-S-MH** tests for **gene-wide episodic positive selection** on a
predefined set of branches.

Key characteristics:

- Gene-level hypothesis test
- Determines whether **any site on any tested branch experienced positive selection**
- Incorporates **synonymous rate variation and multi-hit substitutions**

BUSTED-type methods are often used as a **first-pass test** to determine whether
a gene contains evidence of episodic adaptation before conducting site-level
analyses.

Documentation:  
https://hyphy.org/methods/selection-methods/#busted

---

### CFEL — Contrast Fixed Effects Likelihood

**CFEL** compares **selection pressures between predefined groups of branches**.

Key characteristics:

- Tests whether **site-specific selection differs between two lineages**
- Identifies **lineage-specific evolutionary constraints or adaptations**
- Useful in comparative evolutionary studies

For example, CFEL can test whether a site experiences stronger purifying
selection in one clade compared to another.

Documentation:  
https://hyphy.org/methods/selection-methods/#cfel

---

### RELAX — Selection Intensity Analysis

**RELAX** tests whether **selection has intensified or relaxed** along specific
lineages.

Key characteristics:

- Quantifies shifts in selection strength using parameter **k**
- **k > 1 → intensified selection**
- **k < 1 → relaxed selection**

RELAX is particularly useful for studying evolutionary scenarios such as:

- host shifts  
- changes in population size  
- functional constraint loss  

Documentation:  
https://hyphy.org/methods/selection-methods/#relax

---

### Why multiple tests?

Different evolutionary processes leave different statistical signatures.
AOC integrates multiple complementary approaches to capture these signals.

| Signal | Method |
|------|------|
| Persistent site-level selection | FEL |
| Episodic site-level adaptation | MEME |
| Adaptive lineages | aBSREL |
| Gene-wide episodic adaptation | BUSTED |
| Lineage-specific site differences | CFEL |
| Selection intensity changes | RELAX |

Together, these analyses provide a **comprehensive evolutionary profile of
protein-coding genes**.

---

## Output Files and Interpretation

AOC produces a structured set of outputs that summarize evolutionary selection analyses performed with HyPhy. While HyPhy generates detailed JSON output files for each method, AOC automatically parses these results into tabular summaries and visualizations that are easier to interpret and use for downstream analysis.

The outputs are organized by **sample** and **partition**, allowing users to examine selection signals across gene partitions or alignment segments.

---

### Output Directory Structure

```
results/
  {sample}/
    selection/
      part1/
        FEL.json
        MEME.json
        ABSREL.json
        BUSTEDS-MH.json
        RELAX.json
        CFEL.json
      part2/
        ...
    tables/
      part1/
        {sample}.part1.AOC.FEL_Results.csv
        {sample}.part1.AOC.MEME_Results.csv
        {sample}.part1.AOC.ABSREL_Results.csv
        {sample}.part1.AOC.BUSTEDS-MH_Results.csv
        {sample}.part1.AOC.RELAX_Results.csv
        {sample}.part1.AOC.CFEL_Results.csv
      {sample}.AOC.merged_FEL_Results.csv
      {sample}.AOC.merged_MEME_Results.csv
      {sample}.AOC.merged_ABSREL_Results.csv
      {sample}.AOC.merged_BUSTEDS-MH_Results.csv
      {sample}.AOC.merged_RELAX_Results.csv
      {sample}.AOC.merged_CFEL_Results.csv
      {sample}.selection_overview.csv
    visualizations/
      FEL.merged.png
      MEME.merged.png
```

---

### HyPhy JSON Output

Each selection method produces a **JSON file** containing the full statistical output from HyPhy. These files include:

- likelihood estimates
- substitution rate parameters
- site or branch level statistics
- likelihood ratio test statistics
- p-values and corrected p-values

These JSON files preserve the complete analysis output and can be used for advanced downstream analysis or reproducibility.

HyPhy documentation describing these outputs can be found here:

https://hyphy.org/methods/selection-methods/

Because these JSON files contain nested data structures, AOC automatically converts them into more user-friendly tables.

---

### Tabular Results

For each partition, AOC generates CSV tables summarizing the key statistics from each selection method.

### FEL Results

File:

```
{sample}.partX.AOC.FEL_Results.csv
```

FEL detects **pervasive selection at individual codon sites**.

Typical columns include:

| Column | Meaning |
|------|------|
| CodonSite | Codon position in the alignment |
| alpha | Synonymous substitution rate |
| beta | Nonsynonymous substitution rate |
| dN/dS | Ratio of nonsynonymous to synonymous substitutions |
| p-value | Significance test for selection |
| adjusted_p-value | Multiple testing corrected p-value |

Interpretation:

- **dN/dS > 1** suggests positive selection  
- **dN/dS < 1** suggests purifying selection  
- Sites with **adjusted p-value ≤ 0.10** are typically considered significant.

### MEME Results

File:

```
{sample}.partX.AOC.MEME_Results.csv
```

MEME detects **episodic positive selection** at individual codon sites.

Important columns:

| Column | Meaning |
|------|------|
| CodonSite | Codon position |
| alpha | Synonymous substitution rate |
| beta+ | Nonsynonymous rate under selection |
| p-value | Test for episodic selection |

Interpretation:

- Significant p-values indicate **sites experiencing positive selection on at least one branch of the phylogeny**.

### aBSREL Results

File:

```
{sample}.partX.AOC.ABSREL_Results.csv
```

aBSREL identifies **branches of the phylogeny experiencing episodic diversification**.

Columns include:

| Column | Meaning |
|------|------|
| Branch | Branch name in the phylogenetic tree |
| Corrected P-value | Multiple-testing corrected significance |
| omega_max | Maximum estimated dN/dS rate |
| significant_branch_0.10 | Indicator for branches under selection |

Interpretation:

- Branches with **Corrected P-value ≤ 0.10** show evidence of episodic adaptive evolution.

### BUSTED-S-MH Results

File:

```
{sample}.partX.AOC.BUSTEDS-MH_Results.csv
```

BUSTED-S-MH tests for **gene-wide episodic positive selection**.

Columns include:

| Column | Meaning |
|------|------|
| p_value | Significance of gene-wide selection |
| LRT | Likelihood ratio test statistic |
| tested_branches | Number of foreground branches tested |

Interpretation:

- **p_value ≤ 0.05–0.10** indicates evidence that **at least one site on at least one tested branch experienced positive selection**.

### RELAX Results

File:

```
{sample}.partX.AOC.RELAX_Results.csv
```

RELAX tests for **changes in the intensity of natural selection**.

Columns include:

| Column | Meaning |
|------|------|
| k | Selection intensity parameter |
| p_value | Significance test |
| selection_shift | Relaxed or intensified selection |

Interpretation:

- **k > 1** → intensified selection  
- **k < 1** → relaxed selection  

### CFEL Results

File:

```
{sample}.partX.AOC.CFEL_Results.csv
```

CFEL compares **selection pressures between two groups of branches**.

Columns include:

| Column | Meaning |
|------|------|
| CodonSite | Codon position |
| p-value | Statistical test for differential selection |
| significant_site_0.10 | Indicator for significant differences |

Interpretation:

- Significant sites indicate **different evolutionary pressures between the compared branch sets**.

## Merged Results

For each method, AOC combines partition-level tables into a **single merged file**:

```
{sample}.AOC.merged_FEL_Results.csv
{sample}.AOC.merged_MEME_Results.csv
{sample}.AOC.merged_ABSREL_Results.csv
{sample}.AOC.merged_BUSTEDS-MH_Results.csv
{sample}.AOC.merged_RELAX_Results.csv
{sample}.AOC.merged_CFEL_Results.csv
```

These files include a **Partition column** so results can be compared across partitions.

---

## Selection Overview Table

File:

```
{sample}.selection_overview.csv
```

This table provides a concise summary of the selection signals detected across all analyses.

Example:

| sample | partition | method | metric | value |
|------|------|------|------|------|
| BDNF | 1 | FEL | significant_sites_FDR_0.10 | 3 |
| BDNF | 1 | MEME | significant_sites_FDR_0.10 | 1 |
| BDNF | 1 | ABSREL | significant_branches | 2 |
| BDNF | 1 | RELAX | k | 0.85 |

This overview allows users to quickly identify partitions showing strong signals of selection.

---

## Visualizations

AOC also generates plots summarizing selection signals across the alignment.

Examples include:

```
FEL.merged.png
MEME.merged.png
```

These plots typically display:

- codon site position on the x-axis
- significance or selection metrics on the y-axis

They help visually identify clusters of sites under selection.

---

## Recommended Interpretation Workflow

A typical workflow for interpreting AOC results is:

1. Examine **selection_overview.csv** to identify partitions with strong signals.
2. Inspect **merged FEL and MEME tables** to identify specific codon sites under selection.
3. Use **ABSREL results** to determine which phylogenetic branches experienced adaptive evolution.
4. Evaluate **BUSTED-S-MH** results to determine whether the gene shows gene-wide episodic selection.
5. Examine **RELAX results** to detect shifts in selection intensity.
6. Use visualizations to identify patterns of selection across the alignment.

---

## Additional Documentation

Detailed explanations of each method and its statistical output are available in the HyPhy documentation:

https://hyphy.org/methods/selection-methods/

Users interested in advanced interpretation or methodological details should consult the original HyPhy publications associated with each method.

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
