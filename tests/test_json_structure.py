import json
from pathlib import Path

OUTDIR = Path("tests/tmp_output")

def test_selection_json_valid():
    """
    Ensure selection JSON files contain expected keys.
    """
    selection_dir = OUTDIR / "selection" / "tiny"
    json_files = list(selection_dir.glob("*.json"))

    assert len(json_files) > 0, "No selection JSON files produced"

    for jf in json_files:
        with open(jf) as f:
            data = json.load(f)

        # Check minimal required structure
        assert isinstance(data, dict)
        assert "analysis" in data or "fits" in data