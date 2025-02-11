from collections import defaultdict

import pandas as pd
from metrics.aggregations import basic_aggregations
from metrics.groupings import by_key, main_accident_cause
from metrics.metric import Metric
from utils import save_json

category_series_blueprint = {
    "druhy_nehod": Metric([], by_key("p6"), basic_aggregations),
    "druhy_srazek_jedoucich_vozidel": Metric(
        [], by_key("p7", neq=0), basic_aggregations
    ),
    "druhy_pevnych_prekazek": Metric([], by_key("p8", neq=0), basic_aggregations),
    "zavineni_nehody": Metric(
        [],
        by_key("p10"),
        basic_aggregations,
        ensure_categories=[0, 1, 2, 3, 4, 5, 6, 7],
    ),
    "alkohol_u_vinika": Metric(
        [],
        by_key("p11"),
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
        by_key("p12"),
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
    "druhy_zvere": Metric(
        [],
        by_key("p8a", neq=0),
        basic_aggregations,
        ensure_categories=[
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
        ],
    ),
    "drogy_u_vinika": Metric(
        [],
        by_key("p11a"),
        basic_aggregations,
        ensure_categories=[0, 1, 2, 3, 4, 5, 6, 7, 8],
    ),
}

calculated_series: dict[
    str,
    dict[str, dict[int, list[tuple[int, int, int | float]]]],
] = {key: defaultdict(lambda: defaultdict(list)) for key in category_series_blueprint}

accidents = pd.read_feather(snakemake.input["accidents"])
vehicles = pd.read_feather(snakemake.input["vehicles"])
effects = pd.read_feather(snakemake.input["effects"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])
gps = pd.read_feather(snakemake.input["gps"])

accidents["datetime"] = pd.to_datetime(accidents["p2a"], format="%d.%m.%Y")

for series_name, series_metric in category_series_blueprint.items():
    result = series_metric.apply(
        accidents=accidents,
        vehicles=vehicles,
        effects=effects,
        pedestrians=pedestrians,
        gps=gps,
    )
    for aggregation, data in result.items():
        for (year, month, category), value in data.items():
            calculated_series[series_name][aggregation][int(category)].append(
                (year, month, value)
            )

save_json(snakemake.output[0], calculated_series)
