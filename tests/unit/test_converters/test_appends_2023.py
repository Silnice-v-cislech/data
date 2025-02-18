import numpy as np
import pandas as pd

from workflow.scripts.metrics.converters.appends_2023 import append_first_vehicle


def test_append_first_vehicle():
    accidents_mock = pd.DataFrame(
        data={
            "p1": [0, 1, 2],
        }
    )
    vehicles_mock = pd.DataFrame(
        data={
            "p1": [0, 0, 1],
            "id_vozidla": [1, 2, 1],
            "p44": [3, 8, 9],
        }
    )

    result = append_first_vehicle(accidents=accidents_mock, vehicles=vehicles_mock)

    assert np.array_equal(
        result["accidents"]["p44"], [3.0, 9.0, np.nan], equal_nan=True
    )
