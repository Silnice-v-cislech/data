from typing import TypeAlias
from numbers import Number
from collections import ChainMap
from dataclasses import dataclass
from collections.abc import Callable
import pandas as pd

filter: TypeAlias = Callable[[pd.DataFrame], dict[str, pd.DataFrame]]
grouping: TypeAlias = Callable[[pd.DataFrame], pd.api.typing.DataFrameGroupBy]
aggregation: TypeAlias = Callable[
    [pd.api.typing.DataFrameGroupBy], list[tuple[tuple[int, int], Number]]
]


@dataclass
class Metric:
    """Processing steps for a metric."""

    filters: list[filter]
    grouping: grouping
    aggregations: dict[str, aggregation]

    def apply(self, **kwargs):
        arguments = kwargs
        for filter in self.filters:
            result = filter(**arguments)
            arguments = ChainMap(result, kwargs)
        grouped = self.grouping(**arguments)
        return {
            key: aggregation(grouped) for key, aggregation in self.aggregations.items()
        }
