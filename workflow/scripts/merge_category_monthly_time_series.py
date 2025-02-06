from collections import defaultdict

from utils import load_json, save_json

result: dict[str, dict[str, dict[int, list[tuple[int, int, int | float]]]]] = (
    defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
)

for input_file in snakemake.input:
    content = load_json(input_file)

    for name, data in content.items():
        for aggregation, values in data.items():
            for category, series in values.items():
                result[name][aggregation][category].extend(series)

save_json(snakemake.output[0], result)
