from collections import defaultdict

import pandas as pd
from category_monthly_time_series_common_metrics import (
    accident_cause,
    accident_kind_counts,
    culprit_alcohol,
    moving_vehicle_crash_kind_counts,
    solid_obstacle_kind_counts,
    main_accident_cause,
    main_accident_cause_detailed,
)
from utils import save_json

category_series_blueprint = {
    "pocty_druhu_nehod": accident_kind_counts,
    "pocty_druhu_srazek_jedoucich_vozidel": moving_vehicle_crash_kind_counts,
    "pocty_druhu_pevnych_prekazek": solid_obstacle_kind_counts,
    "zavineni_nehody": accident_cause,
    "alkohol_u_vinika": culprit_alcohol,
    "hlavni_priciny_nehody": main_accident_cause,
    "hlavni_priciny_nehody_podrobne": main_accident_cause_detailed,
}

calculated_series: dict[
    str,
    dict[int, list[tuple[int, int, int | float]]],
] = {key: defaultdict(list) for key in category_series_blueprint}

accidents = pd.read_feather(snakemake.input["accidents"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])

accidents["datetime"] = pd.to_datetime(accidents["p2a"], format="%Y-%m-%d")

for series_name, series_function in category_series_blueprint.items():
    result = series_function(
        accidents=accidents,
        pedestrians=pedestrians,
    )
    for (year, month, category), value in result.items():
        calculated_series[series_name][int(category)].append((year, month, value))

save_json(snakemake.output[0], calculated_series)
