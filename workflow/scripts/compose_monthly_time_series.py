from typing import Any

from utils import load_json, load_text_file, save_pretty_json


def format_records(series):
    return {f"{year}-{month:02}": value for year, month, value in series}


def format_records_selected_months(series):
    return {f"{year}": value for year, _, value in series}


NORMALIZED_SUFFIX_NAME = "_normalized"

SELECTED_MONTHS_SUFFIX_NAME = "_selected_months"

base = load_json(snakemake.input["base"])
normalized = load_json(snakemake.input["normalized"])
selected_months = load_json(snakemake.input["selected_months"])
labels = load_json(snakemake.input["labels"])
categorized = load_json(snakemake.input["categorized"])
source = load_text_file(snakemake.input["source"])

result: dict[str, dict[str, dict[str, Any]]] = {}

for name in base:
    result[name] = {
        "data": {},
        "labels": {},
    }
    base_label = labels.get(name, "")

    result[name]["data"]["base"] = format_records(base[name])
    result[name]["labels"]["base"] = base_label

    result[name]["data"]["normalized"] = format_records(normalized[name])
    result[name]["labels"]["normalized"] = base_label + labels.get(
        "normalized_suffix", ""
    )

    for month in range(1, 12 + 1):
        if str(month) in selected_months[name]:
            records = selected_months[name][str(month)]
        else:
            records = []

        key = f"month{int(month):02}"
        result[name]["data"][key] = format_records_selected_months(records)
        result[name]["labels"][key] = base_label + labels.get(f"{key}_suffix", "")

for name, categories in categorized.items():
    result[name] = {
        "data": {},
        "labels": {},
    }
    for category, series in categories.items():
        result[name]["data"][category] = format_records(series)
        result[name]["labels"][category] = labels.get(f"{name}_{category}", "")


result["meta"] = {
    "source": source,
}

save_pretty_json(snakemake.output[0], result)
