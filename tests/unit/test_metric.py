from datetime import datetime

import pandas as pd

from workflow.scripts.metrics.aggregations import (
    accident_count,
    basic_aggregations,
    death_count,
)
from workflow.scripts.metrics.groupings import all, by_key
from workflow.scripts.metrics.metric import Metric


def test_metric():
    metric = Metric([], all, basic_aggregations)

    result = metric.apply(
        accidents=pd.DataFrame(
            data={
                "p1": [0, 1],
                "p13a": [3, 1],
                "p13b": [6, 4],
                "p13c": [5, 7],
                "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
            }
        )
    )

    assert {key: dict(value) for key, value in result.items()} == {
        "nehody": {(2016, 1): 2},
        "umrti": {(2016, 1): 4},
        "tezke_zraneni": {(2016, 1): 10},
        "lehke_zraneni": {(2016, 1): 12},
    }


def test_ensure_categories():
    metric = Metric(
        [],
        by_key("p6"),
        {"nehody": accident_count, "umrti": death_count},
        ensure_categories=[1, 2, 3, 4],
    )

    result = metric.apply(
        accidents=pd.DataFrame(
            data={
                "p1": [0, 1],
                "p13a": [3, 1],
                "p13b": [6, 4],
                "p13c": [5, 7],
                "p6": [1, 3],
                "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
            }
        )
    )

    assert {key: dict(value) for key, value in result.items()} == {
        "nehody": {
            (2016, 1, 1): 1,
            (2016, 1, 2): 0,
            (2016, 1, 3): 1,
            (2016, 1, 4): 0,
        },
        "umrti": {
            (2016, 1, 1): 3,
            (2016, 1, 2): 0,
            (2016, 1, 3): 1,
            (2016, 1, 4): 0,
        },
    }
