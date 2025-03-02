import pandas as pd


def append_first_vehicle(*, accidents: pd.DataFrame, vehicles: pd.DataFrame, **_):
    """Vehicle available for 2016 format.
    In case the accident is caused by a driver, it is the vehicle of the drive,
    otherwise it is a vehicle that was part of the accident.
    """
    return {
        "accidents": accidents.merge(
            vehicles[vehicles["id_vozidla"] == 1], how="left", on="p1"
        )
    }
