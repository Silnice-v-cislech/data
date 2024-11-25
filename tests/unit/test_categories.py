from datetime import datetime

import pandas as pd

from workflow.scripts.category_monthly_time_series_common_metrics import (
    accident_cause,
    accident_kind_counts,
    culprit_alcohol,
    moving_vehicle_crash_kind_counts,
    solid_obstacle_kind_counts,
)


def test_accident_kind_counts():
    data = pd.DataFrame(
        data={
            "p6": [0, 1, 1, 5, 5],
            "datetime": [
                datetime(2016, 1, 1),
                datetime(2016, 2, 1),
                datetime(2016, 2, 2),
                datetime(2016, 5, 1),
                datetime(2016, 6, 2),
            ],
        }
    )

    result = []
    for (year, month, category), value in accident_kind_counts(accidents=data).items():
        result.append((year, month, category, value))

    assert result == [
        (2016, 1, 0, 1),
        (2016, 2, 1, 2),
        (2016, 5, 5, 1),
        (2016, 6, 5, 1),
    ]


def test_moving_vehicle_crash_kind_counts():
    data = pd.DataFrame(
        data={
            "p7": [0, 1, 1, 5, 5],
            "datetime": [
                datetime(2016, 1, 1),
                datetime(2016, 2, 1),
                datetime(2016, 2, 2),
                datetime(2016, 5, 1),
                datetime(2016, 6, 2),
            ],
        }
    )

    result = []
    for (year, month, category), value in moving_vehicle_crash_kind_counts(
        accidents=data
    ).items():
        result.append((year, month, category, value))

    assert result == [(2016, 2, 1, 2), (2016, 5, 5, 1), (2016, 6, 5, 1)]


def test_solid_obstacle_kind_counts():
    data = pd.DataFrame(
        data={
            "p8": [0, 1, 1, 5, 5],
            "datetime": [
                datetime(2016, 1, 1),
                datetime(2016, 2, 1),
                datetime(2016, 2, 2),
                datetime(2016, 5, 1),
                datetime(2016, 6, 2),
            ],
        }
    )

    result = []
    for (year, month, category), value in solid_obstacle_kind_counts(
        accidents=data
    ).items():
        result.append((year, month, category, value))

    assert result == [(2016, 2, 1, 2), (2016, 5, 5, 1), (2016, 6, 5, 1)]


def test_accident_cause():
    data = pd.DataFrame(
        data={
            "p10": [0, 1, 1, 5, 5],
            "datetime": [
                datetime(2016, 2, 1),
                datetime(2016, 2, 1),
                datetime(2016, 2, 2),
                datetime(2016, 5, 1),
                datetime(2016, 5, 2),
            ],
        }
    )

    result = []
    for (year, month, category), value in accident_cause(accidents=data).items():
        result.append((year, month, category, value))

    assert result == [
        (2016, 2, 0, 1),
        (2016, 2, 1, 2),
        (2016, 2, 2, 0),
        (2016, 2, 3, 0),
        (2016, 2, 4, 0),
        (2016, 2, 5, 0),
        (2016, 2, 6, 0),
        (2016, 2, 7, 0),
        (2016, 5, 0, 0),
        (2016, 5, 1, 0),
        (2016, 5, 2, 0),
        (2016, 5, 3, 0),
        (2016, 5, 4, 0),
        (2016, 5, 5, 2),
        (2016, 5, 6, 0),
        (2016, 5, 7, 0),
    ]


def test_culprit_alcohol():
    data = pd.DataFrame(
        data={
            "p11": [0, 1, 1, 5, 5],
            "datetime": [
                datetime(2016, 2, 1),
                datetime(2016, 2, 1),
                datetime(2016, 2, 2),
                datetime(2016, 5, 1),
                datetime(2016, 5, 2),
            ],
        }
    )

    result = []
    for (year, month, category), value in culprit_alcohol(accidents=data).items():
        result.append((year, month, category, value))

    assert result == [
        (2016, 2, 0, 1),
        (2016, 2, 1, 2),
        (2016, 2, 2, 0),
        (2016, 2, 3, 0),
        (2016, 2, 4, 0),
        (2016, 2, 5, 0),
        (2016, 2, 6, 0),
        (2016, 2, 7, 0),
        (2016, 2, 8, 0),
        (2016, 2, 9, 0),
        (2016, 5, 0, 0),
        (2016, 5, 1, 0),
        (2016, 5, 2, 0),
        (2016, 5, 3, 0),
        (2016, 5, 4, 0),
        (2016, 5, 5, 2),
        (2016, 5, 6, 0),
        (2016, 5, 7, 0),
        (2016, 5, 8, 0),
        (2016, 5, 9, 0),
    ]
