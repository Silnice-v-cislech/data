from utils import load_json, save_pretty_json

def format_records(series):
    return {f"{year}-{month:02}": value for year, month, value in series}

def format_records_selected_months(series):
    return {f"{year}": value for year, _, value in series}

NORMALIZED_SUFFIX_NAME = "_normalized"

SELECTED_MONTHS_SUFFIX_NAME = "_selected_months"

base = load_json(snakemake.input["base"])
normalized = load_json(snakemake.input["normalized"])
selected_months = load_json(snakemake.input["selected_months"])
styles = load_json(snakemake.input["styles"])
labels = load_json(snakemake.input["labels"])

result = {}

for name in base:
    result[name] = {
        "data": {},
        "styles": {},
        "labels": {},
    }
    base_style = styles.get(name, dict())
    base_label = labels.get(name, "")

    result[name]["data"]["base"] = format_records(base[name])
    result[name]["styles"]["base"] = base_style
    result[name]["labels"]["base"] = base_label

    result[name]["data"]["normalized"] = format_records(normalized[name])
    result[name]["styles"]["normalized"] = base_style | styles.get("normalized_update", dict())
    result[name]["labels"]["normalized"] = base_label + labels.get("normalized_suffix", "")

    for month in range(1, 12 + 1):
        if str(month) in selected_months[name]:
            records = selected_months[name][str(month)]
        else:
            records = []
        
        key = f"month{int(month):02}"
        result[name]["data"][key] = format_records_selected_months(records)
        result[name]["styles"][key] = base_style | styles.get(f"{key}_update", dict())
        result[name]["labels"][key] = base_label + labels.get(f"{key}_suffix", "")

save_pretty_json(snakemake.output[0], result)
