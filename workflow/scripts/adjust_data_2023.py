import pandas as pd

accidents = pd.read_feather(snakemake.input["accidents"])
vehicles = pd.read_feather(snakemake.input["vehicles"])
effects = pd.read_feather(snakemake.input["effects"])
pedestrians = pd.read_feather(snakemake.input["pedestrians"])
gps = pd.read_feather(snakemake.input["gps"])

if not accidents.empty:
    # add proper datetime
    accidents["datetime"] = pd.to_datetime(accidents["p2a"], format="%d.%m.%Y")

    # convert last two digits of year to full year
    vehicles = vehicles.join(accidents[["datetime", "p1"]].set_index("p1"), on="p1")

    def complete_year(row):
        digits = row["p47"]
        year = row["datetime"].year
        try:
            digits = int(digits)
        except (ValueError, TypeError):
            return 0
        if digits <= year % 100:  # probably in 2000s
            return 2000 + digits
        else:  # probably in 1900s
            return 1900 + digits

    vehicles["p47"] = vehicles.apply(complete_year, axis=1)
    vehicles = vehicles.drop(columns=["datetime"])

accidents.to_feather(snakemake.output["accidents"])
vehicles.to_feather(snakemake.output["vehicles"])
effects.to_feather(snakemake.output["effects"])
pedestrians.to_feather(snakemake.output["pedestrians"])
gps.to_feather(snakemake.output["gps"])
