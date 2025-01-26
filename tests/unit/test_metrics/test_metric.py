from datetime import datetime
import pandas as pd
from workflow.scripts.metrics.metric import Metric
from workflow.scripts.metrics.aggregations import basic_aggregations
from workflow.scripts.metrics.groupings import all


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
