from calendar import monthrange

from workflow.scripts.utils import load_json, save_json

data = load_json(snakemake.input[0])

for series in data.values():
    for record in series:
        year, month, value = record

        record[2] = value / monthrange(year, month)[1] * 30

save_json(snakemake.output[0], data)
