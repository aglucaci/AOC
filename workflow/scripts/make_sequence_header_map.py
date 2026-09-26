#!/usr/bin/env python3
"""
make_sequence_header_map_and_test_list.py

UPDATED VERSION (FASTA-free)

This version no longer requires a FASTA file. It simply converts the
sequence names in the labels CSV into the normalized FASTA-safe format
used by the AOC pipeline and outputs:

  1) mapping CSV
  2) Test-set text file (one sequence per line)

Normalization rule:
  - replace non-alphanumeric characters with "_"
  - collapse multiple "_"
  - strip leading/trailing "_"

Inputs:
  --labels-csv   e.g. BDNF.small.sequence_labels.csv
  --out-map      e.g. BDNF.sequence_header_map.csv
  --out-test     e.g. BDNF.test_sequences.txt

Expected labels CSV columns:
  label, fasta_sequence_header
(but the script auto-detects the first two columns if names differ)
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------
def norm_header(s: str) -> str:
    """Convert header into FASTA-safe format."""
    s = s.strip()
    s = re.sub(r"[^A-Za-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s)
    return s.strip("_")


# ---------------------------------------------------------------------
# Column detection
# ---------------------------------------------------------------------
def detect_label_and_header_cols(fieldnames: List[str]) -> Tuple[str, str]:
    cols = [c.strip() for c in fieldnames]
    label_col = "label" if "label" in cols else cols[0]
    hdr_col = "fasta_sequence_header" if "fasta_sequence_header" in cols else cols[1]
    return label_col, hdr_col


# ---------------------------------------------------------------------
# Build mapping
# ---------------------------------------------------------------------
def build_mapping(labels_csv: Path):
    rows_out = []
    mapping: Dict[str, str] = {}

    with labels_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or len(reader.fieldnames) < 2:
            raise ValueError("labels CSV must have at least 2 columns")

        label_col, hdr_col = detect_label_and_header_cols(reader.fieldnames)

        for r in reader:
            lab = (r.get(label_col) or "").strip()
            orig = (r.get(hdr_col) or "").strip()

            normalized = norm_header(orig)

            rows_out.append(
                (lab, orig, normalized, normalized, "ok")
            )

            mapping[orig] = normalized

    return rows_out, mapping


# ---------------------------------------------------------------------
# Write mapping CSV
# ---------------------------------------------------------------------
def write_mapping_csv(rows_out, out_csv: Path):
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "label",
            "original_header",
            "fasta_header",
            "normalized_candidate",
            "status",
        ])
        w.writerows(rows_out)


# ---------------------------------------------------------------------
# Write Test list
# ---------------------------------------------------------------------
def write_test_list(labels_csv: Path, mapping: Dict[str, str], out_txt: Path) -> int:
    test_headers: List[str] = []

    with labels_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        label_col, hdr_col = detect_label_and_header_cols(reader.fieldnames)

        for r in reader:
            if (r.get(label_col) or "").strip().lower() == "test":
                orig = (r.get(hdr_col) or "").strip()
                mapped = mapping.get(orig)
                if mapped:
                    test_headers.append(mapped)

    out_txt.parent.mkdir(parents=True, exist_ok=True)
    with out_txt.open("w") as f:
        for h in test_headers:
            f.write(h + "\n")

    return len(test_headers)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Build normalized header mapping CSV and Test-set list (FASTA-free)."
    )
    ap.add_argument("--labels-csv", required=True)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--out-test", required=True)
    args = ap.parse_args()

    labels_csv = Path(args.labels_csv)
    out_map = Path(args.out_map)
    out_test = Path(args.out_test)

    if not labels_csv.exists():
        raise SystemExit(f"[ERROR] labels CSV not found: {labels_csv}")

    rows_out, mapping = build_mapping(labels_csv)

    write_mapping_csv(rows_out, out_map)
    n_test = write_test_list(labels_csv, mapping, out_test)

    print(f"[OK] Mapping CSV written: {out_map}")
    print(f"[OK] Test list written:   {out_test} (n={n_test})")


if __name__ == "__main__":
    main()
