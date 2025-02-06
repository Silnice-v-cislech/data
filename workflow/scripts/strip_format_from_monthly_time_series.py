from utils import load_json, save_pretty_json

data = load_json(snakemake.input[0])

for name, record in data.items():
    if name != "meta":
        record.pop("labels")

save_pretty_json(snakemake.output[0], data)
