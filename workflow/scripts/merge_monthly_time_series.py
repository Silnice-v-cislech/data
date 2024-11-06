from collections import defaultdict

from utils import load_json, save_json

result = defaultdict(list)

for input_file in snakemake.input:
    data = load_json(input_file)

    for key, value in data.items():
        result[key].extend(value)

save_json(snakemake.output[0], result)
