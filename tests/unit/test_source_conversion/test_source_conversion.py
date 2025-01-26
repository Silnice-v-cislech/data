import os
import shutil
import subprocess as sp
import sys
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))


def test_workflow():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath("tests/unit/test_source_conversion/data")
        expected_path = PurePosixPath("tests/unit/test_source_conversion/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "resources/2016/chodci.feather",
                "resources/2016/nehody.feather",
                "resources/2023/chodci.feather",
                "resources/2023/nehody.feather",
                "resources/2023/gps.feather",
                "resources/2023/nasledky.feather",
                "resources/2023/vozidla.feather",
                "-j1",
                "--target-files-omit-workdir-adjustment",
                "--directory",
                workdir,
            ]
        )

        # Check results.
        for path in [
            "2016/nehody.feather",
            "2016/chodci.feather",
            "2023/nehody.feather",
            "2023/vozidla.feather",
            "2023/nasledky.feather",
            "2023/chodci.feather",
            "2023/gps.feather",
        ]:
            assert (workdir / "resources" / path).is_file()

            result = pd.read_feather(workdir / "resources" / path)
            expected = pd.read_feather(expected_path / path)

            try:
                assert result.compare(expected).empty
            except ValueError as e:
                raise AssertionError from e
