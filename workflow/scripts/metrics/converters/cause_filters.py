import pandas as pd


def caused_by_motor_vehicle_driver(*, accidents: pd.DataFrame, **_):
    return {
        "accidents": accidents[accidents["p10"] == 1]
    }  # ZAVINĚNÍ NEHODY == řidičem motorového vozidla
