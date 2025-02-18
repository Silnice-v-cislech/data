from datetime import datetime

import pandas as pd

accidents = pd.read_feather(snakemake.input["accidents"])
vehicles = pd.read_feather(snakemake.input["vehicles"])
effects = pd.read_feather(snakemake.input["effects"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])
gps = pd.read_feather(snakemake.input["gps"])

# add proper datetime
accidents["datetime"] = pd.to_datetime(accidents["p2a"], format="%d.%m.%Y")


# convert last two digits of year to full year
def complete_year(digits):
    try:
        digits = int(digits)
    except (ValueError, TypeError):
        return 0
    if digits <= datetime.now().year % 100:  # probably in 2000s
        return 2000 + digits
    else:  # probably in 1900s
        return 1900 + digits


vehicles["p47"] = vehicles["p47"].apply(complete_year)

accidents.to_feather(snakemake.output["accidents"])
vehicles.to_feather(snakemake.output["vehicles"])
effects.to_feather(snakemake.output["effects"])
pedestrians.to_feather(snakemake.output["pedestrians"])
gps.to_feather(snakemake.output["gps"])
