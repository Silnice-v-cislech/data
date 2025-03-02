from datetime import datetime
from typing import Any

import orjson


def load_text_file(filename):
    with open(filename, "r") as f:
        return f.read().strip()


def load_json(filename):
    with open(filename, "br") as f:
        return orjson.loads(f.read())


def save_json(filename, data):
    with open(filename, "bw") as f:
        f.write(orjson.dumps(data, option=orjson.OPT_NON_STR_KEYS))


def save_pretty_json(filename, data):
    with open(filename, "bw") as f:
        f.write(
            orjson.dumps(
                data,
                option=orjson.OPT_NON_STR_KEYS
                | orjson.OPT_INDENT_2
                | orjson.OPT_SORT_KEYS,
            )
        )


def zero_monthly_time_series(
    start: datetime, end: datetime, value: Any = 0
) -> dict[tuple[int, int], int]:
    start_year = int(start.year)
    start_month = int(start.month)
    end_year = int(end.year)
    end_month = int(end.month)

    if start_year == end_year:
        return {
            (start_year, month): value for month in range(start_month, end_month + 1)
        }

    return (
        {(start_year, month): value for month in range(start_month, 12 + 1)}
        | {
            (year, month): value
            for year in range(start_year + 1, end_year)
            for month in range(1, 12 + 1)
        }
        | {(end_year, month): value for month in range(1, end_month + 1)}
    )
