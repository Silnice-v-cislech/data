import pandas as pd
from utils import save_json
from monthly_time_series_common import accident_count, death_count, severe_injury_count, light_injury_count

series_blueprint = {
    "pocet_nehod": accident_count,
    "pocet_usmrcenych": death_count,
    "tezce_zranenych_osob": severe_injury_count,
    "lehce_zranenych_osob": light_injury_count,
}

calculated_series = {key: [] for key in series_blueprint}

accidents = pd.read_feather(snakemake.input["accidents"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])

accidents["datetime"] = pd.to_datetime(accidents["p2a"], format="%Y-%m-%d")

for series_name, series_function in series_blueprint.items():
    result = series_function(
        accidents=accidents,
        pedestrians=pedestrians,
    )
    for (year, month), value in result.items():
        calculated_series[series_name].append((year, month, value))

save_json(snakemake.output[0], calculated_series)
