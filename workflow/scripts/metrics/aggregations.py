def accident_count(groupby):
    return groupby["datetime"].count()


def death_count(groupby):
    return groupby["p13a"].sum().astype(int)


def severe_injury_count(groupby):
    return groupby["p13b"].sum().astype(int)


def light_injury_count(groupby):
    return groupby["p13c"].sum().astype(int)


basic_aggregations = {
    "nehody": accident_count,
    "umrti": death_count,
    "tezke_zraneni": severe_injury_count,
    "lehke_zraneni": light_injury_count,
}
