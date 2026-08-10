# Developer and Technical Notes

This document is for contributors who want to understand how AOC is organized,
how to make changes safely, and how the workflow fits together internally.
User-facing installation and execution instructions live in `README.md`; issue
reporting and community expectations live in `CONTRIBUTING.md`.

## Contributing

### Development setup

Start from a clean checkout and install the AOC environment:

```bash
bash install.sh AOC envs/AOC.yaml
conda activate AOC
```

Run the smoke test before making changes so you know the local environment is
healthy:

```bash
bash tests/test_installation.sh
```

For a dry run against the default sample sheet:

```bash
bash run_AOC.sh --samples samples.csv --dry-run
```

### Contribution workflow

1. Open an issue or discussion for substantial workflow, dependency, output
   schema, or scientific-method changes.
2. Create a topic branch from the current development branch.
3. Keep changes scoped. Avoid mixing workflow behavior, documentation updates,
   and broad refactors in one pull request unless they are tightly coupled.
4. Update tests or add a small fixture when changing parsing, path resolution,
   branch-label behavior, summarization logic, or command generation.
5. Run the relevant checks and include the commands in the pull request.
6. Document any changed input columns, config keys, output paths, or method
   behavior in `README.md` and this file as appropriate.

### Coding conventions

- Prefer explicit output paths: shared sequence artifacts under
  `results/sequences/{sequence}/` and sample-specific artifacts under
  `results/{sample}/`.
- Keep Snakemake rules deterministic: declare all important inputs and outputs,
  write logs, and avoid relying on hidden files.
- Use `config/config.yaml` for paths, executable names, and tunable parameters
  that users may need to override.
- Keep helper scripts in `scripts/` when logic becomes easier to test outside
  the Snakefile.
- Preserve backward compatibility for existing `samples.csv` files where
  possible. If compatibility cannot be preserved, update the README and include
  a migration note.
- Keep generated outputs out of version control unless they are small,
  intentional fixtures.

### Tests and validation

The current test suite is intentionally small and centered on installation and
stable parser behavior.

```bash
python -m unittest tests/test_busteds_mh.py
bash tests/test_installation.sh
```

Use `tests/data/tiny_samples.csv` and `tests/data/tiny.fasta` for quick workflow
checks. Larger biological examples belong under `data/` only when they are
small enough to be practical for repository users.

When changing the workflow, validate at least one labeled and one unlabeled
sample if the change touches branch labels, RELAX, CFEL, tree labeling, or
selection-method command arguments.

### Pull request checklist

- The affected workflow path has been run or dry-run.
- New or changed outputs are documented.
- Config defaults remain sensible for local runs and Slurm runs.
- Logs are written for new rules or command wrappers.
- Tests cover changed parser or summarization behavior.
- The pull request explains any scientific assumption, statistical threshold,
  or HyPhy option that changed.

## Software Architecture

### High-level design

AOC is a Snakemake-driven analysis pipeline. The repository keeps the main
workflow in one Snakefile and uses small supporting scripts, HyPhy batch files,
and shell wrappers around external tools.

The core flow is:

```text
samples.csv
  -> codon-aware alignment and cleanup
  -> duplicate removal and TN93 downsampling
  -> GARD recombination analysis
  -> dynamic partition FASTA generation
  -> per-partition tree inference and optional branch labeling
  -> HyPhy selection analyses
  -> tables, plots, merged summaries, and HTML reports
```

### Repository layout

- `workflow/Snakefile`: main workflow definition and most result parsing.
- `config/config.yaml`: default runtime configuration, tool paths, output root,
  GARD settings, TN93 settings, and MACSE memory settings.
- `run_AOC.sh`: local command-line wrapper around Snakemake.
- `submit_AOC.slurm`: example Slurm entry point.
- `install.sh`: environment installer and smoke-test helper.
- `envs/AOC.yaml`: conda environment specification.
- `scripts/`: Python and HyPhy helper scripts used by the workflow.
- `software/`: bundled HyPhy utilities such as LabelTrees and duplicate removal.
- `tests/`: smoke tests, unit tests, and tiny workflow fixtures.
- `data/`: small example datasets.
- `JOSS/`: manuscript and publication assets.

### Configuration and inputs

AOC supports two input modes. The legacy mode keeps the FASTA path directly in
`samples.csv`:

```text
sample,codon_fasta,sequence_labels_csv
```

The shared-sequence mode separates foreground test samples from reusable FASTA
inputs:

```text
sample,sequence,sequence_labels_csv
```

paired with `sequences.csv`:

```text
sequence,codon_fasta
```

The Snakefile reads this file through the `samples_csv` config key. Each row
defines a sample-specific execution branch. In shared-sequence mode, the
`sequence` value identifies the shared preprocessing/GARD/tree branch. The
`sample` value becomes the sample-specific output directory name under the
configured output root, which defaults to `results/`.

If `sequence_labels_csv` is omitted or blank, the sample is treated as
unlabeled. Labeled samples are used to build a header map and a Test-branch list
for branch-aware HyPhy methods.

Important config keys include:

- `samples_csv`: sample sheet path.
- `sequences_csv`: sequence-to-FASTA sheet path for shared-sequence mode.
- `outdir`: output root.
- `hyphy`, `hyphy_mpi`, `mpirun`, `fasttree`, `macse_launcher`: executable
  names or paths.
- `gard_processors`: MPI ranks for GARD.
- `tn93_threshold`, `tn93_max_seqs`: TN93 clustering controls.
- `preserve_test_sequences`: protects labeled Test records during duplicate
  removal and TN93 clustering.
- `require_test_sequences`: whether labeled samples must produce Test branches.
- `macse_*`: Java heap and memory tuning for MACSE.

### Workflow stages

#### 1. Preprocessing

The preprocessing rules produce a cleaned codon alignment for each sequence:

1. `macse`: runs MACSE for codon-aware alignment and amino-acid output.
2. `cln`: runs HyPhy CLN on the codon alignment.
3. `strike_ambigs_msa`: removes or masks ambiguous alignment content.
4. `remove_duplicates_msa`: removes duplicate sequences with the bundled HyPhy
   duplicate-removal utility. When `preserve_test_sequences` is enabled, any
   missing Test records are restored from the pre-deduplicated alignment.
5. `tn93_cluster`: clusters and limits retained sequences before recombination
   analysis while retaining protected Test records.

These steps write sequence-local files under `results/sequences/{sequence}/` and
logs under `results/logs/`.

#### 2. Recombination and dynamic partitions

`recombination` runs HyPhy GARD through MPI and produces both the GARD JSON and
`.best-gard` breakpoint file.

The `parse_gard` checkpoint reads the `.best-gard` file, creates codon-aligned
partition FASTA files, and writes a manifest under:

```text
results/sequences/{sequence}/gard/segments/
results/sequences/{sequence}/gard/segments.json
```

This checkpoint is what makes downstream partition fan-out dynamic. Helper
functions read the sequence manifest to discover the actual partitions for each
sample that points at that sequence.

#### 3. Branch-label preparation

For samples with a labels CSV, `build_header_map_and_test_list` sanitizes labels
and sequence headers into HyPhy-safe names. It writes:

```text
results/{sample}/{sample}.sequence_header_map.csv
results/{sample}/{sample}.test_sequences.txt
```

If a sample has no labels, helper functions create or use an empty test-list
file. RELAX and CFEL then emit skipped JSON placeholders instead of failing.

#### 4. Per-partition phylogenetics and selection analyses

For each GARD partition:

1. `fasttree_partition` builds one shared nucleotide tree per sequence partition
   with FastTree.
2. `label_tree_TestForeground_partition` labels Test branches per sample when a non-empty
   test list exists.
3. Selection rules run HyPhy methods:
   - `FEL_partition`
   - `MEME_partition`
   - `ABSREL_partition`
   - `BUSTEDS_MH_partition`
   - `RELAX_partition`
   - `CFEL_partition`

FEL, MEME, aBSREL, and BUSTED-S-MH can run without labels. RELAX and CFEL
require branch labels and are skipped with explicit placeholder JSON when no
Test list exists.

#### 5. Tables, plots, and summaries

Partition JSON files are converted into CSVs under:

```text
results/{sample}/tables/part{part}/
```

FEL and MEME also produce partition-level plots under:

```text
results/{sample}/visualizations/part{part}/
```

Merge rules combine partition-level tables into sample-level tables:

```text
results/{sample}/tables/{sample}.AOC.merged_*_Results.csv
```

`SelectionOverviewTable` builds a compact cross-method overview, and
`executive_summary` writes:

```text
results/{sample}/summary/executive_summary.html
```

### Output contract

The output layout is part of AOC's practical API. Downstream users may rely on:

- `results/{sample}/selection/part{part}/*.json`
- `results/{sample}/tables/part{part}/*.csv`
- `results/{sample}/tables/{sample}.AOC.merged_*_Results.csv`
- `results/{sample}/tables/{sample}.selection_overview.csv`
- `results/{sample}/visualizations/{sample}.FEL.merged.png`
- `results/{sample}/visualizations/{sample}.MEME.merged.png`
- `results/{sample}/summary/executive_summary.html`
- `results/sequences/{sequence}/gard/segments.json`
- `results/sequences/{sequence}/trees/part{part}/FastTree.treefile`

If a change renames, removes, or changes the schema of any of these files, treat
it as a user-facing change and document it.

### External tool boundary

AOC delegates scientific heavy lifting to established external tools:

- MACSE for codon-aware alignment.
- HyPhy for cleaning, GARD, branch labeling helpers, and selection analyses.
- FastTree for per-partition tree inference.
- TN93-based clustering via `scripts/tn93_cluster.py`.
- Altair and `vl-convert-python` for local plot rendering.

The Snakefile is responsible for orchestration, path conventions, branching
logic, parsing, and result aggregation. Avoid reimplementing external method
behavior in Python unless the goal is to parse, validate, or summarize outputs.

### Common extension points

To add a new selection method:

1. Add a per-partition rule that consumes `_partition_fasta` and the labeled
   tree.
2. Add its JSON output to `_all_partition_jsons_for_sample`.
3. Add a table parser rule if the method needs CSV output.
4. Add merge and overview logic if the result should appear in sample-level
   summaries.
5. Add README documentation for method purpose, output files, and interpretation.
6. Add a tiny parser unit test if the JSON structure is nontrivial.

To add a new config option:

1. Add a default and comment in `config/config.yaml`.
2. Read it through `config.get(...)` in the Snakefile or helper script.
3. Document when users should change it.
4. Keep the default compatible with local runs.

To change input sample metadata:

1. Update the required-column validation near the top of `workflow/Snakefile`.
2. Update `README.md` input examples.
3. Add or update tiny fixtures in `tests/data/`.
4. Validate both labeled and unlabeled sample rows.

### Design constraints

- The workflow should remain runnable from the repository root with
  `bash run_AOC.sh --samples samples.csv`.
- Outputs should keep the sequence/sample split: shared sequence artifacts under
  `results/sequences/{sequence}/`, sample-specific artifacts under
  `results/{sample}/`, and selection outputs partition-aware.
- Dynamic partitions should continue to flow through the GARD manifest rather
  than hard-coded partition counts.
- Labeled and unlabeled samples should both be first-class modes.
- HPC support should avoid assumptions that break local execution.
- Logs should be sufficient to diagnose external-tool failures without rerunning
  with ad hoc debug flags.
