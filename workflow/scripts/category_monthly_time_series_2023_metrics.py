from category_monthly_time_series_common_metrics import ensure_categories


def animal_obstacle_kind_counts(*, accidents, **_):
    return ensure_categories(
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22],
        accidents[accidents["p8a"] != 0]
        .groupby(
            [
                accidents["datetime"].dt.year,
                accidents["datetime"].dt.month,
                accidents["p8a"],
            ]
        )["datetime"]
        .count(),
    )


def culprit_drugs(*, accidents, **_):
    return ensure_categories(
        [0, 1, 2, 3, 4, 5, 6, 7, 8],
        accidents.groupby(
            [
                accidents["datetime"].dt.year,
                accidents["datetime"].dt.month,
                accidents["p11a"],
            ]
        )["datetime"].count(),
    )
