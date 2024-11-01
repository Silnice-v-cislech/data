from workflow.scripts.utils import load_json, save_json
from collections import defaultdict

result = defaultdict(list)

for input_file in snakemake.input:
    data = load_json(input_file)
    
    for key, value in data.items():
        result[key].extend(value)

save_json(snakemake.output[0], result)
