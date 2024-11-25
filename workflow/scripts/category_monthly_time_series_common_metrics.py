import pandas as pd


def ensure_categories(categories, data):
    levels = [
        data.index.levels[0].values,
        data.index.levels[1].values,
        categories,
    ]
    new_index = pd.MultiIndex.from_product(levels, names=data.index.names)
    return data.reindex(new_index, fill_value=0)


def accident_kind_counts(*, accidents, **_):
    return accidents.groupby(
        [accidents["datetime"].dt.year, accidents["datetime"].dt.month, accidents["p6"]]
    )["datetime"].count()


def moving_vehicle_crash_kind_counts(*, accidents, **_):
    return (
        accidents[accidents["p7"] != 0]
        .groupby(
            [
                accidents["datetime"].dt.year,
                accidents["datetime"].dt.month,
                accidents["p7"],
            ]
        )["datetime"]
        .count()
    )


def solid_obstacle_kind_counts(*, accidents, **_):
    return (
        accidents[accidents["p8"] != 0]
        .groupby(
            [
                accidents["datetime"].dt.year,
                accidents["datetime"].dt.month,
                accidents["p8"],
            ]
        )["datetime"]
        .count()
    )


def accident_cause(*, accidents, **_):
    return ensure_categories(
        [0, 1, 2, 3, 4, 5, 6, 7],
        accidents.groupby(
            [
                accidents["datetime"].dt.year,
                accidents["datetime"].dt.month,
                accidents["p10"],
            ]
        )["datetime"].count(),
    )


def culprit_alcohol(*, accidents, **_):
    return ensure_categories(
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        accidents.groupby(
            [
                accidents["datetime"].dt.year,
                accidents["datetime"].dt.month,
                accidents["p11"],
            ]
        )["datetime"].count(),
    )
