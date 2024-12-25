from datetime import datetime

import pandas as pd

from workflow.scripts.category_monthly_time_series_common_metrics import (
    accident_cause,
    accident_kind_counts,
    culprit_alcohol,
    moving_vehicle_crash_kind_counts,
    solid_obstacle_kind_counts,
    main_accident_cause,
    main_accident_cause_detailed,
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


def test_main_accident_cause():
    data = pd.DataFrame(
        data={
            "p12": [201, 205, 209, 501, 100],
            "datetime": [
                datetime(2016, 2, 1),
                datetime(2016, 2, 1),
                datetime(2016, 2, 2),
                datetime(2016, 2, 3),
                datetime(2016, 2, 4),
            ],
        }
    )

    result = []
    for (year, month, category), value in main_accident_cause_detailed(accidents=data).items():
        result.append((year, month, category, value))

    assert result == [
        (2016, 2, 100, 1),
        (2016, 2, 201, 1),
        (2016, 2, 202, 0),
        (2016, 2, 203, 0),
        (2016, 2, 204, 0),
        (2016, 2, 205, 1),
        (2016, 2, 206, 0),
        (2016, 2, 207, 0),
        (2016, 2, 208, 0),
        (2016, 2, 209, 1),
        (2016, 2, 301, 0),
        (2016, 2, 302, 0),
        (2016, 2, 303, 0),
        (2016, 2, 304, 0),
        (2016, 2, 305, 0),
        (2016, 2, 306, 0),
        (2016, 2, 307, 0),
        (2016, 2, 308, 0),
        (2016, 2, 309, 0),
        (2016, 2, 310, 0),
        (2016, 2, 311, 0),
        (2016, 2, 401, 0),
        (2016, 2, 402, 0),
        (2016, 2, 403, 0),
        (2016, 2, 404, 0),
        (2016, 2, 405, 0),
        (2016, 2, 406, 0),
        (2016, 2, 407, 0),
        (2016, 2, 408, 0),
        (2016, 2, 409, 0),
        (2016, 2, 410, 0),
        (2016, 2, 411, 0),
        (2016, 2, 412, 0),
        (2016, 2, 413, 0),
        (2016, 2, 414, 0),
        (2016, 2, 501, 1),
        (2016, 2, 502, 0),
        (2016, 2, 503, 0),
        (2016, 2, 504, 0),
        (2016, 2, 505, 0),
        (2016, 2, 506, 0),
        (2016, 2, 507, 0),
        (2016, 2, 508, 0),
        (2016, 2, 509, 0),
        (2016, 2, 510, 0),
        (2016, 2, 511, 0),
        (2016, 2, 512, 0),
        (2016, 2, 513, 0),
        (2016, 2, 514, 0),
        (2016, 2, 515, 0),
        (2016, 2, 516, 0),
        (2016, 2, 601, 0),
        (2016, 2, 602, 0),
        (2016, 2, 603, 0),
        (2016, 2, 604, 0),
        (2016, 2, 605, 0),
        (2016, 2, 606, 0),
        (2016, 2, 607, 0),
        (2016, 2, 608, 0),
        (2016, 2, 609, 0),
        (2016, 2, 610, 0),
        (2016, 2, 611, 0),
        (2016, 2, 612, 0),
        (2016, 2, 613, 0),
        (2016, 2, 614, 0),
        (2016, 2, 615, 0),
    ]

    result = []
    for (year, month, category), value in main_accident_cause(accidents=data).items():
        result.append((year, month, category, value))

    assert result == [
        (2016, 2, 1, 1),
        (2016, 2, 2, 3),
        (2016, 2, 3, 0),
        (2016, 2, 4, 0),
        (2016, 2, 5, 1),
        (2016, 2, 6, 0),
    ]
