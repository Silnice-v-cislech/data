from datetime import datetime

import pandas as pd
import pytest

from workflow.scripts.metrics.aggregations import (
    accident_count,
    death_count,
    light_injury_count,
    severe_injury_count,
)


@pytest.mark.parametrize(
    ("aggregation", "data", "expected"),
    [
        (
            accident_count,
            {
                "p1": [0, 1],
                "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
            },
            {2016: 2},
        ),
        (
            death_count,
            {
                "p1": [0, 1],
                "p13a": [3, 1],
                "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
            },
            {2016: 4},
        ),
        (
            severe_injury_count,
            {
                "p1": [0, 1],
                "p13b": [6, 4],
                "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
            },
            {2016: 10},
        ),
        (
            light_injury_count,
            {
                "p1": [0, 1],
                "p13c": [5, 7],
                "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
            },
            {2016: 12},
        ),
    ],
)
def test_aggregation(aggregation, data, expected):
    dataframe = pd.DataFrame(data=data)
    grouped = dataframe.groupby([dataframe["datetime"].dt.year])
    result = aggregation(grouped)
    assert dict(result) == expected
