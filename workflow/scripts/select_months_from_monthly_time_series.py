from collections import defaultdict
from workflow.scripts.utils import load_json, save_json

data = load_json(snakemake.input[0])

result: dict[str, dict[int, list[int|float]]] = {}

for name, series in data.items():
    result[name] = defaultdict(list)
    for record in series:
        _, month, _ = record
        result[name][month].append(record)

save_json(snakemake.output[0], result)
