import pandas as pd


def all(*, accidents: pd.DataFrame, **_):
    return accidents.groupby(
        [accidents["datetime"].dt.year, accidents["datetime"].dt.month]
    )


def accident_kind(*, accidents, **_):
    return accidents.groupby(
        [accidents["datetime"].dt.year, accidents["datetime"].dt.month, accidents["p6"]]
    )


def moving_vehicle_crash_kind(*, accidents, **_):
    return accidents[accidents["p7"] != 0].groupby(
        [
            accidents["datetime"].dt.year,
            accidents["datetime"].dt.month,
            accidents["p7"],
        ]
    )


def solid_obstacle_kind(*, accidents, **_):
    return accidents[accidents["p8"] != 0].groupby(
        [
            accidents["datetime"].dt.year,
            accidents["datetime"].dt.month,
            accidents["p8"],
        ]
    )


def accident_cause(*, accidents, **_):
    return accidents.groupby(
        [
            accidents["datetime"].dt.year,
            accidents["datetime"].dt.month,
            accidents["p10"],
        ]
    )


def culprit_alcohol(*, accidents, **_):
    return accidents.groupby(
        [
            accidents["datetime"].dt.year,
            accidents["datetime"].dt.month,
            accidents["p11"],
        ]
    )


def main_accident_cause_detailed(*, accidents, **_):
    return accidents.groupby(
        [
            accidents["datetime"].dt.year,
            accidents["datetime"].dt.month,
            accidents["p12"],
        ]
    )


def main_accident_cause(*, accidents, **_):
    return accidents.groupby(
        [
            accidents["datetime"].dt.year,
            accidents["datetime"].dt.month,
            accidents["p12"].apply(lambda x: x // 100),
        ]
    )


def animal_obstacle_kind(*, accidents, **_):
    return accidents[accidents["p8a"] != 0].groupby(
        [
            accidents["datetime"].dt.year,
            accidents["datetime"].dt.month,
            accidents["p8a"],
        ]
    )


def culprit_drugs(*, accidents, **_):
    return accidents.groupby(
        [
            accidents["datetime"].dt.year,
            accidents["datetime"].dt.month,
            accidents["p11a"],
        ]
    )
