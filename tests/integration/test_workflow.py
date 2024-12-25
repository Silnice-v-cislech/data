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

        EMPTY_2016_FORMAT = {
            "nehody": {
                "p2a": [],
                "p6": [],
                "p7": [],
                "p8": [],
                "p10": [],
                "p11": [],
                "p12": [],
                "p13a": [],
                "p13b": [],
                "p13c": [],
            },
            "chodci": {},
        }
        EMPTY_2023_FORMAT = {
            "nehody": {
                "p2a": [],
                "p6": [],
                "p7": [],
                "p8": [],
                "p8a": [],
                "p10": [],
                "p11": [],
                "p11a": [],
                "p12": [],
                "p13a": [],
                "p13b": [],
                "p13c": [],
            },
            "chodci": {},
            "nasledky": {},
            "vozidla": {},
            "gps": {},
        }
        # create input files
        data_2016_format = {
            2016: {
                "nehody": {
                    "p1": [0, 1, 2],
                    "p2a": ["2016-01-01", "2016-02-01", "2016-02-02"],
                    "p6": [1, 0, 3],
                    "p7": [1, 0, 2],
                    "p8": [1, 0, 2],
                    "p10": [1, 0, 2],
                    "p11": [1, 0, 2],
                    "p12": [100, 202, 209],
                    "p13a": [2, 3, 0],
                    "p13b": [4, 5, 0],
                    "p13c": [10, 15, 0],
                },
                "chodci": {},
            },
            2017: EMPTY_2016_FORMAT,
            2018: EMPTY_2016_FORMAT,
            2019: EMPTY_2016_FORMAT,
            2020: EMPTY_2016_FORMAT,
            2021: EMPTY_2016_FORMAT,
            2022: EMPTY_2016_FORMAT,
        }
        data_2023_format = {
            2023: {
                "nehody": {
                    "p1": [0, 1, 2],
                    "p2a": ["01.01.2023", "02.01.2023", "02.02.2023"],
                    "p6": [2, 2, 2],
                    "p7": [3, 3, 3],
                    "p8": [8, 8, 8],
                    "p8a": [1, 1, 0],
                    "p10": [7, 7, 6],
                    "p11": [8, 8, 9],
                    "p11a": [6, 6, 0],
                    "p12": [501, 516, 100],
                    "p13a": [3, 5, 0],
                    "p13b": [2, 3, 0],
                    "p13c": [12, 13, 0],
                },
                "chodci": {},
                "nasledky": {},
                "vozidla": {},
                "gps": {},
            },
            2024: EMPTY_2023_FORMAT,
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

        with open(workdir / "results" / "formatted_monthly_time_series.json", "r") as f:
            result = f.read()

        with open(expected_path / "formatted_monthly_time_series.json", "r") as f:
            expected = f.read()

        assert result == expected
