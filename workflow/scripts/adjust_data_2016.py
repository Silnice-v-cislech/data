from datetime import datetime

import pandas as pd

accidents = pd.read_feather(snakemake.input["accidents"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])

# add proper datetime
accidents["datetime"] = pd.to_datetime(accidents["p2a"], format="%Y-%m-%d")


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


accidents["p47"] = accidents["p47"].apply(complete_year)

accidents.to_feather(snakemake.output["accidents"])
pedestrians.to_feather(snakemake.output["pedestrians"])
