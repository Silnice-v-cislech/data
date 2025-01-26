import pandas as pd
from datetime import datetime
from workflow.scripts.metrics.groupings import all


def test_all():
    dataframe = pd.DataFrame(
        data={
            "p1": [0, 1],
            "datetime": [datetime(2016, 1, 1), datetime(2016, 1, 1)],
        }
    )

    result = all(accidents=dataframe)

    assert isinstance(result, pd.api.typing.DataFrameGroupBy)

    assert dict(result["datetime"].count()) == {(2016, 1): 2}
