import numpy as np
import pandas as pd

from workflow.scripts.metrics.converters.cause_filters import (
    caused_by_motor_vehicle_driver,
)


def test_append_first_vehicle():
    accidents_mock = pd.DataFrame(
        data={
            "p1": [0, 1, 2],
            "p10": [1, 2, 1],
        }
    )
    result = caused_by_motor_vehicle_driver(accidents=accidents_mock)

    assert np.array_equal(result["accidents"]["p1"], [0, 2])
