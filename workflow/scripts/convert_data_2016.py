import io
import rarfile  # type: ignore
import pandas as pd

archive = snakemake.input[0]

files_nehody = [
    "00.csv",
    "01.csv",
    "02.csv",
    "03.csv",
    "04.csv",
    "05.csv",
    "06.csv",
    "07.csv",
    "08.csv",
    "09.csv",
    "10.csv",
    "11.csv",
    "12.csv",
    "13.csv",
    "14.csv",
    "15.csv",
    "16.csv",
    "17.csv",
    "18.csv",
    "19.csv",
]
file_chodci = "CHODCI.csv"

nehody_codes = [
    "p1",
    "p36",
    "p37",
    "p2a",
    "weekday(p2a)",
    "p2b",
    "p6",
    "p7",
    "p8",
    "p9",
    "p10",
    "p11",
    "p12",
    "p13a",
    "p13b",
    "p13c",
    "p14",
    "p15",
    "p16",
    "p17",
    "p18",
    "p19",
    "p20",
    "p21",
    "p22",
    "p23",
    "p24",
    "p27",
    "p28",
    "p34",
    "p35",
    "p39",
    "p44",
    "p45a",
    "p47",
    "p48a",
    "p49",
    "p50a",
    "p50b",
    "p51",
    "p52",
    "p53",
    "p55a",
    "p57",
    "p58",
    "a",
    "b",
    "d",
    "e",
    "f",
    "g",
    "h",
    "i",
    "j",
    "k",
    "l",
    "n",
    "o",
    "p",
    "q",
    "r",
    "s",
    "t",
    "p5a",
]
nehody_dtypes: dict[str, str | type[str]] = {
    "p1": str,
    "p36": "Int32",
    "p37": "Int32",
    "p2a": str,
    "weekday(p2a)": "Int32",
    "p2b": str,
    "p6": "Int32",
    "p7": "Int32",
    "p8": "Int32",
    "p9": "Int32",
    "p10": "Int32",
    "p11": "Int32",
    "p12": "Int32",
    "p13a": "Int32",
    "p13b": "Int32",
    "p13c": "Int32",
    "p14": "Int32",
    "p15": "Int32",
    "p16": "Int32",
    "p17": "Int32",
    "p18": "Int32",
    "p19": "Int32",
    "p20": "Int32",
    "p21": "Int32",
    "p22": "Int32",
    "p23": "Int32",
    "p24": "Int32",
    "p27": "Int32",
    "p28": "Int32",
    "p34": "Int32",
    "p35": "Int32",
    "p39": "Int32",
    "p44": "Int32",
    "p45a": "Int32",
    "p47": "Int32",
    "p48a": "Int32",
    "p49": "Int32",
    "p50a": "Int32",
    "p50b": "Int32",
    "p51": "Int32",
    "p52": "Int32",
    "p53": "Int32",
    "p55a": "Int32",
    "p57": "Int32",
    "p58": "Int32",
    "a": str,
    "b": str,
    "d": str,
    "e": str,
    "f": str,
    "g": str,
    "h": str,
    "i": str,
    "j": str,
    "k": str,
    "l": str,
    "n": str,
    "o": str,
    "p": str,
    "q": str,
    "r": str,
    "s": str,
    "t": str,
    "p5a": "Int32",
}

chodci_codes = [
    "p1",
    "p29",
    "p30",
    "p31",
    "p32",
]
chodci_dtypes: dict[str, str | type[str]] = {
    "p1": str,
    "p29": "Int32",
    "p30": "Int32",
    "p31": "Int32",
    "p32": "Int32",
}


with rarfile.RarFile(archive) as rf:
    loaded = []
    for file in files_nehody:
        try:
            if not (
                df := pd.read_csv(
                    io.BytesIO(rf.read(file)),
                    sep=";",
                    encoding="cp1250",
                    names=nehody_codes,
                    dtype=nehody_dtypes,
                    na_values=["", "XX", "--", "AN"],
                )
            ).empty:
                loaded.append(df)
        except Exception as e:
            raise type(e)(f"While processing {archive}: {file}") from e

    df = pd.concat(loaded, ignore_index=True)
    df.to_feather(snakemake.output["accidents"])

    df = pd.read_csv(
        io.BytesIO(rf.read(file_chodci)),
        sep=";",
        encoding="cp1250",
        names=chodci_codes,
        dtype=chodci_dtypes,
    )
    df.to_feather(snakemake.output["pedestrians"])
