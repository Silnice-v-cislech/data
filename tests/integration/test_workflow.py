import os
import re
import shutil
import subprocess as sp
import sys
from itertools import chain
from pathlib import Path, PurePosixPath

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))

RESOURCES_DIR = "resources"


def test_workflow():

    workdir = Path("tests") / "integration" / "temp_test_dir"
    shutil.rmtree(workdir, ignore_errors=True)

    resources_path = PurePosixPath("tests/integration/data") / RESOURCES_DIR
    expected_path = PurePosixPath("tests/integration/expected")

    # Copy data to the temporary workdir.
    shutil.copytree(resources_path.parent, workdir)

    EMPTY_2016_FORMAT = {
        "nehody": {
            "p1": [],
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
            "p47": [],
        },
        "chodci": {},
    }
    EMPTY_2023_FORMAT = {
        "nehody": {
            "p1": [],
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
        "vozidla": {
            "p1": [],
            "id_vozidla": [],
            "p47": [],
        },
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
                "p47": [22, 44, np.nan],
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
                "p1": [0, 1, 2, 3],
                "p2a": ["01.01.2023", "02.01.2023", "02.02.2023", "02.03.2023"],
                "p6": [2, 2, 2, 1],
                "p7": [3, 3, 3, 1],
                "p8": [8, 8, 8, 0],
                "p8a": [1, 1, 0, 0],
                "p10": [7, 7, 6, 1],
                "p11": [8, 8, 9, 0],
                "p11a": [6, 6, 0, 0],
                "p12": [501, 516, 100, 201],
                "p13a": [3, 5, 0, 7],
                "p13b": [2, 3, 0, 5],
                "p13c": [12, 13, 0, 17],
            },
            "chodci": {},
            "nasledky": {},
            "vozidla": {
                "p1": [0, 1, 2, 3, 3, 3],
                "id_vozidla": [1, 1, 1, 1, 2, 3],
                "p47": [0, np.nan, None, np.nan, 34, 89],
            },
            "gps": {},
        },
        2024: EMPTY_2023_FORMAT,
        2025: EMPTY_2023_FORMAT,
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

    # remove processed timestamp
    result = re.sub(r"\"processed\": \".*?\"", '"processed": ""', result)

    assert result == expected

    with open(workdir / "results" / "formatted_monthly_time_series.json", "r") as f:
        result = f.read()

    with open(expected_path / "formatted_monthly_time_series.json", "r") as f:
        expected = f.read()

    # remove processed timestamp
    result = re.sub(r"\"processed\": \".*?\"", '"processed": ""', result)

    assert result == expected
