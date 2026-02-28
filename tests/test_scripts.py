from scripts.write_manifest import main as write_manifest_main
from pathlib import Path
import tempfile

def test_manifest_script_runs():
    """
    Ensures manifest script runs without crashing.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        out = Path(tmpdir) / "manifest.csv"
        write_manifest_main(
            samples_csv="tests/data/tiny_samples.csv",
            out=str(out),
            outdir="tests/tmp_output"
        )
        assert out.exists()