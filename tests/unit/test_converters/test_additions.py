from datetime import datetime

import numpy as np
import pandas as pd

from workflow.scripts.metrics.converters.additions import add_vehicle_age_to_accidents


def test_add_vehicle_age_to_accidents():
    accidents_mock = pd.DataFrame(
        data={
            "p1": [0, 1, 2],
            "p47": [2016, 2024, 0],
            "datetime": [
                datetime(2024, 1, 10),
                datetime(2024, 1, 11),
                datetime(2024, 1, 12),
            ],
        }
    )

    result = add_vehicle_age_to_accidents(accidents=accidents_mock)

    assert np.array_equal(
        result["accidents"]["vehicle_age"], [8.0, 0.0, np.nan], equal_nan=True
    )
