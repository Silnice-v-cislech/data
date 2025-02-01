from collections import defaultdict

import pandas as pd
from metrics.aggregations import basic_aggregations
from metrics.groupings import all
from metrics.metric import Metric
from utils import save_json

series_blueprint = {
    "celkovy_pocet": Metric([], all, basic_aggregations),
}

calculated_series: dict[str, dict[str, list[tuple[int, int, int | float]]]] = {
    key: defaultdict(list) for key in series_blueprint
}

accidents = pd.read_feather(snakemake.input["accidents"])
vehicles = pd.read_feather(snakemake.input["vehicles"])
effects = pd.read_feather(snakemake.input["effects"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])
gps = pd.read_feather(snakemake.input["gps"])

accidents["datetime"] = pd.to_datetime(accidents["p2a"], format="%d.%m.%Y")

for series_name, series_metric in series_blueprint.items():
    result = series_metric.apply(
        accidents=accidents,
        vehicles=vehicles,
        effects=effects,
        pedestrians=pedestrians,
        gps=gps,
    )
    for aggregation, data in result.items():
        for (year, month), value in data.items():
            calculated_series[series_name][aggregation].append((year, month, value))

save_json(snakemake.output[0], calculated_series)
