## Community & Contribution Guidelines

AOC is an actively developed, modular Snakemake-based framework for ortholog discovery, evolutionary analysis, and selection summarization. We welcome community involvement, feature suggestions, and pull requests.

---

## Reporting Issues

Please report bugs or request features via GitHub Issues:

https://github.com/aglucaci/Analysis-of-Orthologous-Collections/issues

When reporting, please include:

- **AOC version or commit hash**
- **Operating system**
- **Conda environment details** (`conda list` or `conda env export`)
- **Command used** (e.g., `snakemake --cores 8 --config samples_csv=...`)
- Relevant logs (especially from `logs/`)
- A minimal reproducible example (e.g., 1 ortholog + `samples.csv`)

If the issue involves selection models (FEL, MEME, aBSREL, BUSTED-S-MH, etc.), please indicate:

- Whether recombination (GARD) was enabled  
- The exact input FASTA and tree (if possible)

---

## Seeking Support

For usage questions or clarification:

- Open a GitHub Issue and label it **question**
- Or contact the maintainer (see Contact section)

Before posting, please verify:

- `pytest` passes locally  
- `bash scripts/test_installation.sh` succeeds  
- Your `samples.csv` is correctly formatted  

---

## Automated Testing

AOC includes both **unit tests and lightweight integration tests**.

### Run the full test suite:

```bash
pip install pytest
pytest