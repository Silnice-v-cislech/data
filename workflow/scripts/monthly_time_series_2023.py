import pandas as pd
from workflow.scripts.utils import save_json
from workflow.scripts.monthly_time_series_common import accident_count, death_count, severe_injury_count, light_injury_count

series_blueprint = {
    "pocet_nehod": accident_count,
    "pocet_usmrcenych": death_count,
    "tezce_zranenych_osob": severe_injury_count,
    "lehce_zranenych_osob": light_injury_count,
}

calculated_series: dict[str, list[tuple[int, int, int|float]]] = {key: [] for key in series_blueprint}

accidents = pd.read_feather(snakemake.input["accidents"])
vehicles = pd.read_feather(snakemake.input["vehicles"])
effects = pd.read_feather(snakemake.input["effects"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])
gps = pd.read_feather(snakemake.input["gps"])

accidents["datetime"] = pd.to_datetime(accidents["p2a"], format="%d.%m.%Y")

for series_name, series_function in series_blueprint.items():
    result = series_function(
        accidents=accidents,
        vehicles=vehicles,
        effects=effects,
        pedestrians=pedestrians,
        gps=gps,
    )
    for (year, month), value in result.items():
        calculated_series[series_name].append((year, month, value))

save_json(snakemake.output[0], calculated_series)
