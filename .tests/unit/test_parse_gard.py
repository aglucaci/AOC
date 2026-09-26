import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

import pytest

sys.path.insert(0, os.path.dirname(__file__))

import common

@pytest.mark.skip(reason="Claims there's no rules to make targets?")
def test_parse_gard():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/config")
        data_path = PurePosixPath(".tests/unit/parse_gard/data")
        expected_path = PurePosixPath(".tests/unit/parse_gard/expected")

        # Copy config to the temporary workdir
        shutil.copytree(config_path, workdir / "config")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir, dirs_exist_ok=True)

        # dbg
        print("results/sequences/BDNF/gard/segments.json results/sequences/BDNF/gard/segments results/sequences/BDNF/BDNF.gard.log", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/sequences/BDNF/gard/segments.json results/sequences/BDNF/gard/segments results/sequences/BDNF/BDNF.gard.log",
            "-f", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
