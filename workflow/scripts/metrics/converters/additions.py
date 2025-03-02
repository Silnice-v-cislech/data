import pandas as pd


def add_vehicle_age_to_accidents(*, accidents: pd.DataFrame, **_):
    accidents["vehicle_age"] = accidents["datetime"].dt.year - accidents["p47"]
    accidents.loc[accidents["p47"] == 0, "vehicle_age"] = None
    return {"accidents": accidents}
