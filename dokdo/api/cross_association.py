import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
import statsmodels.stats.multitest as multi
from qiime2 import Artifact, Metadata

def cross_association_table(
    artifact, target, method='spearman', normalize=None, alpha=0.05,
    multitest='fdr_bh'
):
    """
    Parameters
    ----------
    artifact : str, qiime2.Artifact, or pandas.DataFrame
        Artifact file or object with the semantic type
        ``FeatureTable[Frequency]``. Alternatively, a
        :class:`pandas.DataFrame` object.
    target : pandas.DataFrame
        Target :class:`pandas.DataFrame` object to be used for
        cross-correlation analysis.
    method : {'spearman', 'pearson'}, default: 'spearman'
        Association method.
    normalize : {None, 'log10', 'clr', 'zscore'}, default: None
        Whether to normalize the feature table:

        - None: Do not normalize.
        - 'log10': Apply the log10 transformation after adding a psuedocount
          of 1.
        - 'clr': Apply the centre log ratio (CLR) transformation adding a
          psuedocount of 1.
        - 'zscore': Apply the Z score transformation.
    alpha : float, default: 0.05
        FWER, family-wise error rate.
    multitest : str, default: 'fdr_bh'
        Method used for testing and adjustment of p values, as defined in
        :meth:`statsmodels.stats.multitest.multipletests`.

    Returns
    -------
    pandas.DataFrame
        Cross-association table.
    """
    if method not in ['pearson', 'spearman']:
        raise ValueError("Method must be 'pearson' or 'spearman'")

    funcs = {'pearson': pearsonr, 'spearman': spearmanr}

    if isinstance(artifact, Artifact):
        input_df = artifact.view(pd.DataFrame)
    elif isinstance(artifact, str):
        input_df = Artifact.load(artifact).view(pd.DataFrame)
    elif isinstance(artifact, pd.DataFrame):
        input_df = artifact.copy()
    else:
        raise TypeError(f'Incorrect input type: {type(artifact)}.')

    if normalize == 'log10':
        input_df = input_df.applymap(lambda x: np.log10(x + 1))
    elif normalize == 'clr':
        input_df = input_df.apply(lambda x: clr(x + 1), axis=1, result_type='broadcast')
    elif normalize == 'zscore':
        input_df = input_df.apply(zscore, axis=1, result_type='broadcast')
    elif not normalize:
        pass
    else:
        raise ValueError(f"Not supported normalization: '{normalize}'")

    def compute_corr(x, y):
        results = funcs[method](x, y)
        return results.correlation

    def compute_pval(x, y):
        results = funcs[method](x, y)
        return results.pvalue

    corr_df = input_df.apply(target.corrwith, method=compute_corr)
    pval_df = input_df.apply(target.corrwith, method=compute_pval)

    df1 = corr_df.melt(var_name='taxon', value_name='correlation', ignore_index=False)
    df2 = pval_df.melt(var_name='taxon', value_name='pval', ignore_index=False)
    df3 = pd.concat([df1, df2['pval']], axis=1)
    df3.index.name = 'target'
    df3 = df3.reset_index()
    df3['pvaladj'] = multi.multipletests(df3['pval'], method=multitest)[1]
    df3 = df3.sort_values('pval')

    return df3

def cross_association_heatmap(
    artifact, target, method='spearman', normalize=None, figsize=None,
    cmap=None
):
    """
    Parameters
    ----------
    artifact : str, qiime2.Artifact, pandas.DataFrame
        Artifact file or object with the semantic type
        ``FeatureTable[Frequency]``. Alternatively, a
        :class:`pandas.DataFrame` object.
    """
    df = cross_association_table(
        artifact, target, method=method, normalize=normalize
    )

    df = df.set_index(['target', 'taxon'])['correlation'].unstack()

    g = sns.clustermap(df, figsize=figsize, cmap=cmap)

    return g
