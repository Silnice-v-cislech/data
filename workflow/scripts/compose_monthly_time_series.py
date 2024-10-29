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
format = load_json(snakemake.input["format"])

result = {}

for name, records in base.items():
    result[name] = format.get(name, dict())
    result[name].update({
        "data": {
            "s0": format_records(records)
        }
    })
    
for name, records in normalized.items():
    name_normalized = f"{name}{NORMALIZED_SUFFIX_NAME}"
    result[name_normalized] = format.get(name_normalized, dict())
    result[name_normalized].update({
        "data": {
            "s0": format_records(records)
        }
    })
    
for name, records in selected_months.items():
    name_selected_months = f"{name}{SELECTED_MONTHS_SUFFIX_NAME}"
    result[name_selected_months] = format.get(name_selected_months, dict())
    result[name_selected_months].update({
        "data": {
            f"s{int(key):02}": format_records_selected_months(value)
            for key, value in records.items()
        }
    })
    

save_pretty_json(snakemake.output[0], result)
