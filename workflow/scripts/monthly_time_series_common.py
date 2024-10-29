def accident_count(*, accidents, **_):
    return accidents.groupby(
        [accidents["datetime"].dt.year, accidents["datetime"].dt.month]
    )["datetime"].count()


def death_count(*, accidents, **_):
    return (
        accidents.groupby(
            [accidents["datetime"].dt.year, accidents["datetime"].dt.month]
        )["p13a"]
        .sum()
        .astype(int)
    )


def severe_injury_count(*, accidents, **_):
    return (
        accidents.groupby(
            [accidents["datetime"].dt.year, accidents["datetime"].dt.month]
        )["p13b"]
        .sum()
        .astype(int)
    )


def light_injury_count(*, accidents, **_):
    return (
        accidents.groupby(
            [accidents["datetime"].dt.year, accidents["datetime"].dt.month]
        )["p13c"]
        .sum()
        .astype(int)
    )
