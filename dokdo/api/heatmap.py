from . import common

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from skbio.stats.composition import clr
from qiime2 import Artifact

def heatmap(
    artifact, metadata=None, hue1=None, hue_order1=None,
    hue1_cmap='tab10', hue1_loc='upper right', hue2=None,
    hue_order2=None, hue2_cmap='Pastel1', hue2_loc='upper left',
    normalize=None, method='average', metric='euclidean',
    figsize=(10, 10), row_cluster=True, col_cluster=True, **kwargs
):
    """
    Create a hierarchically clustered heatmap of a feature table.

    Parameters
    ----------
    artifact : str, qiime2.Artifact, pandas.DataFrame
        Artifact file or object with the semantic type
        ``FeatureTable[Frequency]``. Alternatively, a
        :class:`pandas.DataFrame` object.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    hue1 : str, optional
        First grouping variable that will produce labels with different
        colors.
    hue_order1 : list, optional
        Specify the order of categorical levels of the 'hue1' semantic.
    hue1_cmap : str, default: 'tab10'
        Name of the colormap passed to :meth:`matplotlib.cm.get_cmap()` for
        `hue1`.
    hue1_loc : str, default: 'upper right'
        Location of the legend for `hue1`.
    hue2 : str, optional
        Second grouping variable that will produce labels with different
        colors.
    hue_order2 : list, optional
        Specify the order of categorical levels of the 'hue2' semantic.
    hue2_cmap : str, default: 'Pastel1'
        Name of the colormap passed to :meth:`matplotlib.cm.get_cmap()` for
        `hue2`.
    hue2_loc : str, default: 'upper left'
        Location of the legend for `hue2`.
    normalize : {None, 'log10', 'clr'}, default: None
        Normalize the feature table by adding a psuedocount of 1 and then
        taking the log10 of the table or performing centre log ratio
        transformation.
    method : str, default: 'average'
        Linkage method to use for calculating clusters. See
        :meth:`scipy.cluster.hierarchy.linkage()` for more details.
    metric : str, default: 'euclidean'
        Distance metric to use for the data. See
        `scipy.spatial.distance.pdist()` documentation for more options.
    figsize : tuple, default: (10, 10)
        Width, height in inches. Format: (float, float).
    row_cluster : bool, default: True
        If True, cluster the rows.
    col_cluster : bool, default: True
        If True, cluster the columns.
    kwargs : other keyword arguments
        Other keyword arguments will be passed down to
        :meth:`seaborn.clustermap()`.

    Returns
    -------
    seaborn.matrix.ClusterGrid
        A ClusterGrid instance.

    Examples
    --------
    Below is a simple example:

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()
        qza_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/table.qza'
        dokdo.heatmap(qza_file,
                      normalize='log10')

    .. image:: images/heatmap-1.png

    We can color the samples by ``body-site`` and use the centered log-ratio
    transformation (CLR) for normalziation:

    .. code:: python3

        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'
        dokdo.heatmap(qza_file,
                      metadata=metadata_file,
                      normalize='clr',
                      hue1='body-site')

    .. image:: images/heatmap-2.png

    We can omit the clustering of samples:

    .. code:: python3

        dokdo.heatmap(qza_file,
                      metadata=metadata_file,
                      normalize='clr',
                      hue1='body-site',
                      row_cluster=False)

    .. image:: images/heatmap-3.png

    We can add an additional grouping variable ``subject``. Note that
    ``xticklabels`` and ``yticklabels`` are extra keyword arguments that
    are passed to :meth:`seaborn.clustermap`.

    .. code:: python3

        dokdo.heatmap(qza_file,
                      metadata=metadata_file,
                      normalize='clr',
                      hue1='body-site',
                      hue2='subject',
                      xticklabels=False,
                      yticklabels=False)

    .. image:: images/heatmap-4.png

    Finally, we can provide :class:`pandas.DataFrame` as input:

    .. code:: python3

        import pandas as pd
        csv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/table.csv'
        df = pd.read_csv(csv_file, index_col=0)
        dokdo.heatmap(df,
                      metadata=metadata_file,
                      normalize='clr',
                      hue1='body-site')

    .. image:: images/heatmap-5.png
    """
    # Check the input type.
    if isinstance(artifact, Artifact):
        df = artifact.view(pd.DataFrame)
    elif isinstance(artifact, str):
        df = Artifact.load(artifact).view(pd.DataFrame)
    elif isinstance(artifact, pd.DataFrame):
        df = artifact
    else:
        raise TypeError(f'Incorrect feature table type: {type(artifact)}')

    # If the metadata is provided, filter the samples accordingly.
    if metadata is not None:
        mf = common.get_mf(metadata)
        df = pd.concat([df, mf], axis=1, join='inner')
        df.drop(mf.columns, axis=1, inplace=True)
        df = df.loc[:, (df != 0).any(axis=0)]

    # If the hue argument(s) are used, get the row colors.
    lut1 = None
    lut2 = None
    row_colors = None
    if hue1 is not None:
        colors1 = plt.cm.get_cmap(hue1_cmap).colors
        df = pd.concat([df, mf], axis=1, join='inner')
        if hue_order1 is None:
            keys1 = df[hue1].unique()
        else:
            keys1 = hue_order1
            df = df[df[hue1].isin(hue_order1)]
        lut1 = dict(zip(keys1, colors1[:len(keys1)]))
        row_colors = df[hue1].map(lut1)
        df.drop(mf.columns, axis=1, inplace=True)
    if hue2 is not None:
        colors2 = plt.cm.get_cmap(hue2_cmap).colors
        df = pd.concat([df, mf], axis=1, join='inner')
        if hue_order2 is None:
            keys2 = df[hue2].unique()
        else:
            keys2 = hue_order2
            df = df[df[hue2].isin(hue_order2)]
        lut2 = dict(zip(keys2, colors2[:len(keys2)]))
        s = df[hue2].map(lut2)
        row_colors = pd.concat([row_colors, s], axis=1)
        df.drop(mf.columns, axis=1, inplace=True)

    # Apply the appropriate normalziation.
    if normalize == 'log10':
        df = df.apply(lambda x: np.log10(x + 1))
    elif normalize == 'clr':
        df = df.apply(lambda x: clr(x + 1), axis=1, result_type='broadcast')
    else:
        pass

    # Draw the heatmap.
    g = sns.clustermap(df, method=method, metric=metric, figsize=figsize,
                       row_cluster=row_cluster, col_cluster=col_cluster,
                       row_colors=row_colors, **kwargs)

    # If the hue argument(s) are used, add the legend(s).
    if hue1 is not None:
        handles = [Patch(facecolor=lut1[name]) for name in lut1]
        legend1 = plt.legend(handles, lut1, title=hue1, bbox_to_anchor=(1, 1),
                   bbox_transform=plt.gcf().transFigure, loc=hue1_loc)
    if hue2 is not None:
        if hue1 is None:
            raise ValueError("Argument 'hue2' was used without 'hue1'. "
                             "Use 'hue1' instead.")
        handles = [Patch(facecolor=lut2[name]) for name in lut2]
        plt.legend(handles, lut2, title=hue2, bbox_to_anchor=(1, 1),
                   bbox_transform=plt.gcf().transFigure, loc=hue2_loc)
        plt.gca().add_artist(legend1)

    return g
