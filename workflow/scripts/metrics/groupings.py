import pandas as pd


def all(*, accidents: pd.DataFrame, **_):
    return accidents.groupby(
        [accidents["datetime"].dt.year, accidents["datetime"].dt.month]
    )
