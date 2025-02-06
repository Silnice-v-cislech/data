from typing import Any

import pandas as pd


def all(*, accidents: pd.DataFrame, **_):
    return accidents.groupby(
        [accidents["datetime"].dt.year, accidents["datetime"].dt.month]
    )


def by_key(key: str, neq: Any = None):
    if neq is not None:

        def grouping_neq(*, accidents, **_):
            return accidents[accidents[key] != neq].groupby(
                [
                    accidents["datetime"].dt.year,
                    accidents["datetime"].dt.month,
                    accidents[key],
                ]
            )

        return grouping_neq

    def grouping(*, accidents, **_):
        return accidents.groupby(
            [
                accidents["datetime"].dt.year,
                accidents["datetime"].dt.month,
                accidents[key],
            ]
        )

    return grouping


def main_accident_cause(*, accidents, **_):
    return accidents.groupby(
        [
            accidents["datetime"].dt.year,
            accidents["datetime"].dt.month,
            accidents["p12"].apply(lambda x: x // 100),
        ]
    )
