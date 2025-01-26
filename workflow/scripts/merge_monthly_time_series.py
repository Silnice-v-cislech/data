from collections import defaultdict

from utils import load_json, save_json

result = defaultdict(lambda: defaultdict(list))

for input_file in snakemake.input:
    content = load_json(input_file)

    for name, data in content.items():
        for aggregation, values in data.items():
            result[name][aggregation].extend(values)

save_json(snakemake.output[0], result)
