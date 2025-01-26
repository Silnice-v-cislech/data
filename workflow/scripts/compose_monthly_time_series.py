from typing import Any

from utils import load_json, load_text_file, save_pretty_json

from datetime import datetime


def format_records(series):
    return {f"{year}-{month:02}": value for year, month, value in series}


base = load_json(snakemake.input["base"])
labels = load_json(snakemake.input["labels"])
categorized = load_json(snakemake.input["categorized"])
source = load_text_file(snakemake.input["source"])

result: dict[str, dict[str, dict[str, Any]]] = {}

for name in base:
    result[name] = {
        "data": {},
        "labels": {},
    }
    for aggregation, series in base[name].items():
        result[name]["data"][aggregation] = format_records(series)
        result[name]["labels"][aggregation] = labels.get(name, {}).get(aggregation, "")

for name, categories in categorized.items():
    result[name] = {
        "data": {},
        "labels": {},
    }
    for category, series in categories.items():
        result[name]["data"][category] = format_records(series)
        result[name]["labels"][category] = labels.get(f"{name}_{category}", "")


result["meta"] = {"source": source, "processed": datetime.now()}

save_pretty_json(snakemake.output[0], result)
