from collections import defaultdict

from utils import load_json, save_json

result: dict[str, dict[int, list[tuple[int, int, int | float]]]] = defaultdict(
    lambda: defaultdict(list)
)

for input_file in snakemake.input:
    data = load_json(input_file)

    for name, categorized_series in data.items():
        for category, series in categorized_series.items():
            result[name][category].extend(series)

save_json(snakemake.output[0], result)
