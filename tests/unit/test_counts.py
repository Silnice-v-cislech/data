from datetime import datetime

import pandas as pd

from workflow.scripts.monthly_time_series_common import (
    accident_count,
    death_count,
    light_injury_count,
    severe_injury_count,
)

data = {
    "p1": [0, 1],
    "p2a": ["2016-01-01", "2016-02-01"],
    "p13a": [2, 3],
    "p13b": [4, 5],
    "p13c": [10, 15],
}


def test_accident_count():
    data = pd.DataFrame(
        data={
            "p1": [0, 1],
            "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
        }
    )

    result = []
    for (year, month), value in accident_count(accidents=data).items():
        result.append((year, month, value))

    assert result == [(2016, 1, 2)]


def test_death_count():
    data = pd.DataFrame(
        data={
            "p1": [0, 1],
            "p13a": [3, 1],
            "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
        }
    )

    result = []
    for (year, month), value in death_count(accidents=data).items():
        result.append((year, month, value))

    assert result == [(2016, 1, 4)]


def test_severe_injury_count():
    data = pd.DataFrame(
        data={
            "p1": [0, 1],
            "p13b": [6, 4],
            "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
        }
    )

    result = []
    for (year, month), value in severe_injury_count(accidents=data).items():
        result.append((year, month, value))

    assert result == [(2016, 1, 10)]


def test_light_injury_count():
    data = pd.DataFrame(
        data={
            "p1": [0, 1],
            "p13c": [5, 7],
            "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
        }
    )

    result = []
    for (year, month), value in light_injury_count(accidents=data).items():
        result.append((year, month, value))

    assert result == [(2016, 1, 12)]
