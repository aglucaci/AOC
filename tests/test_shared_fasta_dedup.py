import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_byte_identical_fastas_share_preprocessing_jobs(tmp_path):
    snakemake = shutil.which("snakemake")
    if snakemake is None:
        pytest.skip("snakemake is not installed")

    original = REPO_ROOT / "tests" / "data" / "tiny.fasta"
    duplicate = tmp_path / "same-content-different-name.fasta"
    shutil.copyfile(original, duplicate)

    samples_csv = tmp_path / "samples.csv"
    samples_csv.write_text(
        "sample,codon_fasta,sequence_labels_csv\n"
        f"red,{original},\n"
        f"blue,{duplicate},\n",
        encoding="utf-8",
    )

    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT)
    command = [
        snakemake,
        "--snakefile",
        str(REPO_ROOT / "workflow" / "Snakefile"),
        "--directory",
        str(REPO_ROOT),
        "--dry-run",
        "--cores",
        "1",
        "--runtime-source-cache-path",
        str(tmp_path / "source-cache"),
        "--config",
        f"samples_csv={samples_csv}",
        f"outdir={tmp_path / 'results'}",
    ]
    result = subprocess.run(command, env=env, text=True, capture_output=True, check=False)

    assert result.returncode == 0, result.stdout + result.stderr
    output = result.stdout + result.stderr
    for rule in [
        "macse",
        "cln",
        "strike_ambigs_msa",
        "remove_duplicates_msa",
        "tn93_cluster",
        "recombination",
        "parse_gard",
    ]:
        assert re.search(rf"^{re.escape(rule)}\s+1$", output, re.MULTILINE), output

    assert re.search(r"^run_selection_all_partitions\s+2$", output, re.MULTILINE), output
