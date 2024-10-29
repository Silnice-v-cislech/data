import orjson

def load_json(filename):
    with open(filename, "br") as f:
        return orjson.loads(f.read())

def save_json(filename, data):
    with open(filename, "bw") as f:
        f.write(orjson.dumps(data, option=orjson.OPT_NON_STR_KEYS))

def save_pretty_json(filename, data):
    with open(filename, "bw") as f:
        f.write(orjson.dumps(data, option=orjson.OPT_NON_STR_KEYS|orjson.OPT_INDENT_2))