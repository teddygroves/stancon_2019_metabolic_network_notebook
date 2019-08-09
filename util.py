import pandas as pd


def expand_series_of_dicts(s: pd.Series) -> pd.DataFrame:
    """get a dataframe from a series whose values are dictionaries"""
    return pd.DataFrame.from_records(s.values, index=s.index)


def join_string_values(df: pd.DataFrame) -> pd.Series:
    return df.apply(lambda x: '_'.join(x.dropna().astype(str)), axis=1)


def match_string_to_file(s: str, path: str) -> bool:
    if os.path.exists(path):
        return open(path, 'r').read() == s
    else:
        return False

