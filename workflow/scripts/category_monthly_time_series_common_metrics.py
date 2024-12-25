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


def main_accident_cause_detailed(*, accidents, **_):
    return ensure_categories(
        [
            100,
            201,
            202,
            203,
            204,
            205,
            206,
            207,
            208,
            209,
            301,
            302,
            303,
            304,
            305,
            306,
            307,
            308,
            309,
            310,
            311,
            401,
            402,
            403,
            404,
            405,
            406,
            407,
            408,
            409,
            410,
            411,
            412,
            413,
            414,
            501,
            502,
            503,
            504,
            505,
            506,
            507,
            508,
            509,
            510,
            511,
            512,
            513,
            514,
            515,
            516,
            601,
            602,
            603,
            604,
            605,
            606,
            607,
            608,
            609,
            610,
            611,
            612,
            613,
            614,
            615,
        ],
        accidents.groupby(
            [
                accidents["datetime"].dt.year,
                accidents["datetime"].dt.month,
                accidents["p12"],
            ]
        )["datetime"].count(),
    )


def main_accident_cause(*, accidents, **_):
    return ensure_categories(
        [1, 2, 3, 4, 5, 6],
        accidents.groupby(
            [
                accidents["datetime"].dt.year,
                accidents["datetime"].dt.month,
                accidents["p12"].apply(lambda x: x // 100),
            ]
        )["datetime"].count(),
    )
