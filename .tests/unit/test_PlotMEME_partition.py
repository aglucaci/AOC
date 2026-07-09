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
def test_PlotMEME_partition():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/config")
        data_path = PurePosixPath(".tests/unit/PlotMEME_partition/data")
        expected_path = PurePosixPath(".tests/unit/PlotMEME_partition/expected")

        # Copy config to the temporary workdir
        shutil.copytree(config_path, workdir / "config")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir, dirs_exist_ok=True)

        # dbg
        print("results/mini/tables/part1/mini.part1.AOC.MEME_Results.csv results/mini/visualizations/part1/mini.part1.MEME.png results/mini/visualizations/part1/mini.part1.MEME.svg results/mini/visualizations/part1/mini.part1.MEME.figure_legend.txt", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/mini/tables/part1/mini.part1.AOC.MEME_Results.csv results/mini/visualizations/part1/mini.part1.MEME.png results/mini/visualizations/part1/mini.part1.MEME.svg results/mini/visualizations/part1/mini.part1.MEME.figure_legend.txt",
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
