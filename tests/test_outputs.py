import subprocess
from pathlib import Path

OUTDIR = Path("tests/tmp_output")

def test_pipeline_runs():
    """
    Runs minimal ortholog workflow.
    Should complete quickly.
    """
    cmd = [
        "snakemake",
        "--cores", "1",
        "--snakefile", "Snakefile",
        "--config",
        "samples_csv=tests/data/tiny_samples.csv",
        f"outdir={OUTDIR}",
        "--quiet"
    ]
    subprocess.run(cmd, check=True)


def test_expected_files_exist():
    """
    Verify key output files exist.
    """
    expected = [
        OUTDIR / "alignments" / "tiny.fasta",
        OUTDIR / "trees" / "tiny.tree.nwk",
        OUTDIR / "summary" / "run_manifest.csv"
    ]

    for f in expected:
        assert f.exists(), f"Missing expected output file: {f}"