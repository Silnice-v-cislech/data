from collections import defaultdict

import pandas as pd
from metrics.aggregations import basic_aggregations
from metrics.groupings import all
from metrics.metric import Metric
from utils import save_json, zero_monthly_time_series

series_blueprint = {
    "celkovy_pocet": Metric([], all, basic_aggregations),
}

calculated_series: dict[str, dict[str, list[tuple[int, int, int | float]]]] = {
    key: defaultdict(list) for key in series_blueprint
}

accidents = pd.read_feather(snakemake.input["accidents"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])

if not accidents.empty:
    min_date = accidents["datetime"].min()
    max_date = accidents["datetime"].max()

    for series_name, series_metric in series_blueprint.items():
        result = series_metric.apply(
            accidents=accidents,
            pedestrians=pedestrians,
        )
        for aggregation, data in result.items():
            data = zero_monthly_time_series(min_date, max_date) | dict(data.items())
            for (year, month), value in data.items():
                calculated_series[series_name][aggregation].append((year, month, value))

save_json(snakemake.output[0], calculated_series)
