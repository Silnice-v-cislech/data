import os
import shutil
import subprocess as sp
import sys
from itertools import chain
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))

RESOURCES_DIR = "resources"


def test_workflow():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        resources_path = PurePosixPath("tests/integration/data") / RESOURCES_DIR
        expected_path = PurePosixPath("tests/integration/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(resources_path.parent, workdir)

        # create input files
        data_2016_format = {
            2016: {
                "nehody": {
                    "p1": [0, 1],
                    "p2a": ["2016-01-01", "2016-02-01"],
                    "p13a": [2, 3],
                    "p13b": [4, 5],
                    "p13c": [10, 15],
                },
                "chodci": {},
            },
            2017: {
                "nehody": {"p2a": [], "p13a": [], "p13b": [], "p13c": []},
                "chodci": {},
            },
            2018: {
                "nehody": {"p2a": [], "p13a": [], "p13b": [], "p13c": []},
                "chodci": {},
            },
            2019: {
                "nehody": {"p2a": [], "p13a": [], "p13b": [], "p13c": []},
                "chodci": {},
            },
            2020: {
                "nehody": {"p2a": [], "p13a": [], "p13b": [], "p13c": []},
                "chodci": {},
            },
            2021: {
                "nehody": {"p2a": [], "p13a": [], "p13b": [], "p13c": []},
                "chodci": {},
            },
            2022: {
                "nehody": {"p2a": [], "p13a": [], "p13b": [], "p13c": []},
                "chodci": {},
            },
        }
        data_2023_format = {
            2023: {
                "nehody": {
                    "p1": [0, 1],
                    "p2a": ["01.01.2023", "02.01.2023"],
                    "p13a": [3, 5],
                    "p13b": [2, 3],
                    "p13c": [12, 13],
                },
                "chodci": {},
                "nasledky": {},
                "vozidla": {},
                "gps": {},
            },
            2024: {
                "nehody": {"p2a": [], "p13a": [], "p13b": [], "p13c": []},
                "chodci": {},
                "nasledky": {},
                "vozidla": {},
                "gps": {},
            },
        }

        for year, files in chain(data_2016_format.items(), data_2023_format.items()):
            for filename, contents in files.items():
                path = workdir / RESOURCES_DIR / str(year) / (filename + ".feather")
                os.makedirs(path.parent, exist_ok=True)
                df = pd.DataFrame(data=contents)
                df.to_feather(path)

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "results/monthly_time_series.json",
                "-j1",
                "--target-files-omit-workdir-adjustment",
                "--directory",
                workdir,
            ]
        )

        # Check results.
        with open(workdir / "results" / "monthly_time_series.json", "r") as f:
            result = f.read()

        with open(expected_path / "monthly_time_series.json", "r") as f:
            expected = f.read()

        assert result == expected
