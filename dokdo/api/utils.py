import numpy as np
import pandas as pd
from skbio.stats.composition import clr
from qiime2 import Artifact

def import_feature_table(artifact):
    """
    Import given feature table.
    """
    if isinstance(artifact, Artifact):
        df = artifact.view(pd.DataFrame)
    elif isinstance(artifact, str):
        df = Artifact.load(artifact).view(pd.DataFrame)
    elif isinstance(artifact, pd.DataFrame):
        df = artifact
    else:
        raise TypeError(f"Incorrect input type: {type(artifact)}")
    return df

def normalize_feature_table(df, method):
    """
    Normalize given feature table.
    """
    if method == 'log10':
        df = df.applymap(lambda x: np.log10(x + 1))
    elif method == 'clr':
        df = df.apply(lambda x: clr(x + 1), axis=1, result_type='broadcast')
    elif method == 'zscore':
        df = df.apply(zscore, axis=1, result_type='broadcast')
    else:
        raise ValueError(f"Incorrect normalization method: {method}")
    return df

def sort_by_mean(df):
    """
    Sort given feature table by mean taxa abundance.
    """
    _ = df.div(df.sum(axis=1), axis=0)
    _ = _.loc[:, _.mean().sort_values(ascending=False).index]
    return df[_.columns]
