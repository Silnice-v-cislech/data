from collections import ChainMap
from collections.abc import Callable
from dataclasses import dataclass
from numbers import Number
from typing import Any, Optional, TypeAlias

import pandas as pd

converter: TypeAlias = Callable[[pd.DataFrame], dict[str, pd.DataFrame]]
grouping: TypeAlias = Callable[[pd.DataFrame], Any]
aggregation: TypeAlias = Callable[[Any], list[tuple[tuple[int, int], Number]]]


@dataclass
class Metric:
    """Processing steps for a metric."""

    converters: list[converter]
    grouping: grouping
    aggregations: dict[str, aggregation]
    ensure_categories: Optional[list[Any]] = None

    def _apply_ensure_categories(self, data):
        levels = [
            data.index.levels[0].values,
            data.index.levels[1].values,
            self.ensure_categories,
        ]
        new_index = pd.MultiIndex.from_product(levels, names=data.index.names)
        return data.reindex(new_index, fill_value=0)

    def apply(self, **kwargs):
        arguments = kwargs
        for converter in self.converters:
            result = converter(**arguments)
            arguments = ChainMap(result, kwargs)
        grouped = self.grouping(**arguments)
        if self.ensure_categories is None:
            return {
                key: aggregation(grouped)
                for key, aggregation in self.aggregations.items()
            }
        else:
            return {
                key: self._apply_ensure_categories(aggregation(grouped))
                for key, aggregation in self.aggregations.items()
            }
