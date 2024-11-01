import io

import pandas as pd
import rarfile  # type: ignore

archive = snakemake.input[0]

files = {
    "Inehody.xls": snakemake.output["accidents"],
    "IVozidla.xls": snakemake.output["vehicles"],
    "Inasledky.xls": snakemake.output["effects"],
    "Ichodci.xls": snakemake.output["pedestrians"],
    "IntGPS.xls": snakemake.output["gps"],
}

with rarfile.RarFile(archive) as rf:

    for input_file, output_file in files.items():
        try:
            df = pd.read_html(io.BytesIO(rf.read(input_file)))[0]
        except Exception as e:
            raise type(e)(f"While processing {archive}: {input_file}") from e

        df = df.iloc[
            :, :-1
        ]  # "Unnamed" column is added at the end with all NaN, this removes it

        df.to_feather(output_file)
