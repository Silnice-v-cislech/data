from collections import defaultdict

import pandas as pd
from metrics.aggregations import accident_count, basic_aggregations
from metrics.groupings import (
    accident_cause,
    accident_kind,
    culprit_alcohol,
    main_accident_cause,
    main_accident_cause_detailed,
    moving_vehicle_crash_kind,
    solid_obstacle_kind,
)
from metrics.metric import Metric
from utils import save_json

category_series_blueprint = {
    "druhy_nehod": Metric([], accident_kind, basic_aggregations),
    "druhy_srazek_jedoucich_vozidel": Metric(
        [], moving_vehicle_crash_kind, basic_aggregations
    ),
    "druhy_pevnych_prekazek": Metric([], solid_obstacle_kind, basic_aggregations),
    "zavineni_nehody": Metric(
        [],
        accident_cause,
        basic_aggregations,
        ensure_categories=[0, 1, 2, 3, 4, 5, 6, 7],
    ),
    "alkohol_u_vinika": Metric(
        [],
        culprit_alcohol,
        basic_aggregations,
        ensure_categories=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    ),
    "hlavni_priciny_nehody": Metric(
        [],
        main_accident_cause,
        basic_aggregations,
        ensure_categories=[1, 2, 3, 4, 5, 6],
    ),
    "hlavni_priciny_nehody_podrobne": Metric(
        [],
        main_accident_cause_detailed,
        basic_aggregations,
        ensure_categories=[
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
    ),
}

calculated_series: dict[
    str,
    dict[str, dict[int, list[tuple[int, int, int | float]]]],
] = {key: defaultdict(lambda: defaultdict(list)) for key in category_series_blueprint}

accidents = pd.read_feather(snakemake.input["accidents"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])

accidents["datetime"] = pd.to_datetime(accidents["p2a"], format="%Y-%m-%d")

for series_name, series_metric in category_series_blueprint.items():
    result = series_metric.apply(
        accidents=accidents,
        pedestrians=pedestrians,
    )
    for aggregation, data in result.items():
        for (year, month, category), value in data.items():
            calculated_series[series_name][aggregation][int(category)].append(
                (year, month, value)
            )

save_json(snakemake.output[0], calculated_series)
