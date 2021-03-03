from qiime2 import Artifact
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ..api import get_mf
from skbio.stats.composition import clr
from matplotlib.patches import Patch

def heatmap(table, metadata=None, hue=None, hue_order=None, normalize=None,
            method='average', metric='euclidean', figsize=(10, 10),
            row_cluster=True, col_cluster=True, cmap_name='tab10', **kwargs):
    """
    This method creates a heatmap representation of a feature table.

    Parameters
    ----------
    table : str or qiime2.Artifact
        Artifact file or object corresponding to FeatureTable[Frequency].
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    hue : str, optional
        Grouping variable that will produce labels with different colors.
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    normalize : str, optional
        Normalize the feature table by adding a psuedocount of 1 and then
        taking the log10 of the table or performing centre log ratio
        transformation. Choices: {'log10', 'clr'}.
    method : str, default: 'average'
        Linkage method to use for calculating clusters. See
        `scipy.cluster.hierarchy.linkage()` documentation for more
        information.
    metric : str, default: 'euclidean'
        Distance metric to use for the data. See
        `scipy.spatial.distance.pdist()` documentation for more options.
    figsize : tuple, default: (10, 10)
        Width, height in inches. Format: (float, float).
    row_cluster : bool, default: True
        If True, cluster the rows.
    col_cluster : bool, default: True
        If True, cluster the columns.
    cmap_name : str, default: 'tab10'
        Name of the colormap passed to `matplotlib.cm.get_cmap()`.
    kwargs : other keyword arguments
        All other keyword arguments are passed to `seaborn.clustermap`.

    Returns
    -------
    seaborn.matrix.ClusterGrid
        A ClusterGrid instance.
    """
    if isinstance(table, Artifact):
        table = table
    elif isinstance(table, str):
        table = Artifact.load(table)
    else:
        raise TypeError(f"Incorrect feature table type: {type(table)}")

    df = table.view(pd.DataFrame)

    if metadata is not None and hue is not None:
        colors = plt.cm.get_cmap(cmap_name).colors
        mf = get_mf(metadata)
        df = pd.concat([df, mf], axis=1, join='inner')
        df = df.loc[:, (df != 0).any(axis=0)]
        if hue_order is None:
            keys = df[hue].unique()
        else:
            keys = hue_order
            df = df[df[hue].isin(hue_order)]
        lut = dict(zip(keys, colors[:len(keys)]))
        row_colors = df[hue].map(lut)
        df.drop(mf.columns, axis=1, inplace=True)
    elif metadata is None and hue is not None:
        raise ValueError("Argument 'hue' requires 'metadata' argument")
    else:
        lut = None
        row_colors = None

    if normalize == 'log10':
        df = df.apply(lambda x: np.log10(x + 1))
    elif normalize == 'clr':
        df = df.apply(lambda x: clr(x + 1), axis=1, result_type='broadcast')
    else:
        pass

    g = sns.clustermap(df,
                       method=method,
                       metric=metric,
                       figsize=figsize,
                       row_cluster=row_cluster,
                       col_cluster=col_cluster,
                       row_colors=row_colors,
                       **kwargs)

    if lut is not None:
        handles = [Patch(facecolor=lut[name]) for name in lut]
        plt.legend(handles, lut, title=hue,
                   bbox_to_anchor=(1, 1),
                   bbox_transform=plt.gcf().transFigure, loc='upper right')

    return g
