from . import common, utils

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from qiime2 import Artifact, Metadata
from scipy.stats import zscore

def _intersect_samples(df, metadata):
    if metadata is None:
        mf = None
    else:
        mf = common.get_mf(metadata)
        df = pd.concat([df, mf], axis=1, join='inner')
        df = df.drop(mf.columns, axis=1)
        df = df.loc[:, (df != 0).any(axis=0)]
    return df, mf

def heatmap(
    artifact, metadata=None, where=None, sort_samples=None,
    pretty_taxa=False, pname_kws=None, normalize=None, samples=None,
    taxa=None, flip=False, vmin=None, vmax=None, cbar=True, cbar_kws=None,
    cbar_ax=None, square=False, label_columns=None, count=0,
    xticklabels=True, yticklabels=True, ax=None, figsize=None, **kwargs
):
    """
    Create a heatmap representation of a feature table.

    It is strongly recommended to normalize the feature table before
    generating a heatmap with the ``normalize`` option.

    Parameters
    ----------
    artifact : str, qiime2.Artifact, pandas.DataFrame
        Artifact file or object with the semantic type
        ``FeatureTable[Frequency]``. Alternatively, a
        :class:`pandas.DataFrame` object.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    where : str, optional
        SQLite WHERE clause specifying sample metadata criteria that must
        be met to be included in the filtered feature table.
    sort_samples : bool, default: False
        If True, sort the samples by name.
    pretty_taxa : bool, default: False
        If True, display only the smallest taxa rank.
    pname_kws : dict, optional
        Keyword arguments for :meth:`dokdo.api.pname` when ``pretty_taxa``
        is True.
    normalize : {None, 'log10', 'clr', 'zscore'}, default: None
        Whether to normalize the feature table:

        - None: Do not normalize.
        - 'log10': Apply the log10 transformation after adding a psuedocount
          of 1.
        - 'clr': Apply the centre log ratio (CLR) transformation adding a
          psuedocount of 1.
        - 'zscore': Apply the Z score transformation.

    samples, taxa : list, optional
        Specify samples and taxa to be displayed.
    flip : bool, default: False
        If True, flip the x and y axes.
    vmin, vmax : floats, optional
        Values to anchor the colormap, otherwise they are inferred from the
        data and other keyword arguments.
    cbar : bool, default: True
        Whether to draw a colorbar.
    cbar_kws : dict, optional
        Keyword arguments for matplotlib.figure.Figure.colorbar().
    cbar_ax : matplotlib.axes.Axes, optional
        Axes in which to draw the colorbar, otherwise take space from the
        main Axes.
    square : bool, default: False
        If True, set the Axes aspect to 'equal' so each cell will be
        square-shaped.
    label_columns : list, optional
        List of metadata columns to be concatenated to form new sample
        labels. Use the string 'sample-id' to indicate the sample ID column.
    count : int, default: 0
        Number of top taxa to display. When 0, display all.
    xticklabels : bool, default: True
        Whether to plot tick labels for the x-axis.
    yticklabels : bool, default: True
        Whether to plot tick labels for the y-axis.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    dokdo.api.clustermap.clustermap

    Examples
    --------
    Below is a simple example:

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()
        qza_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/table-l3.qza'
        dokdo.heatmap(qza_file,
                      normalize='log10',
                      flip=True,
                      figsize=(10, 7))
        plt.tight_layout()

    .. image:: images/heatmap-1.png

    We can display a heatmap for each sample group:

    .. code:: python3

        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'

        fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(1, 5,
                                                      figsize=(20, 8),
                                                      gridspec_kw={'width_ratios': [1, 1, 1, 1, 0.1]})

        kwargs = dict(normalize='log10',
                      flip=True,
                      linewidths=0.5,
                      metadata=metadata_file,
                      xticklabels=True)

        qza_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/table-l3.qza'

        dokdo.heatmap(qza_file, ax=ax1, where="[body-site] IN ('gut')", cbar=False, yticklabels=True, **kwargs)
        dokdo.heatmap(qza_file, ax=ax2, where="[body-site] IN ('left palm')", yticklabels=False, cbar=False, **kwargs)
        dokdo.heatmap(qza_file, ax=ax3, where="[body-site] IN ('right palm')", yticklabels=False, cbar=False, **kwargs)
        dokdo.heatmap(qza_file, ax=ax4, where="[body-site] IN ('tongue')", yticklabels=False, cbar_ax=ax5, **kwargs)

        ax1.set_title('Gut')
        ax2.set_title('Left palm')
        ax3.set_title('Right palm')
        ax4.set_title('Toungue')

        plt.tight_layout()

    .. image:: images/heatmap-2.png
    """
    df = utils.import_feature_table(artifact)
    df, mf = _intersect_samples(df, metadata)
    if normalize is not None:
        df = utils.normalize_feature_table(df, normalize)

    # Determine which samples to display.
    if where is None and samples is not None:
        df = df.loc[samples]
    elif where is not None and samples is None:
        if isinstance(metadata, str):
            metadata = Metadata.load(metadata)
        samples = metadata.get_ids(where)
        df = df.loc[samples]
    else:
        pass

    # Sort the samples by name, if necessary.
    if sort_samples:
        df = df.sort_index()

    # Determine which taxa to display.
    if taxa is not None and count:
        raise ValueError("Cannot use 'taxa' and 'count' together.")
    elif taxa is not None and not count:
        df = df[taxa]
    elif taxa is None and count:
        cols = df.mean().sort_values(ascending=False).index[:count]
        df = df[cols]
    else:
        pass

    # Flip the axes.
    if flip:
        df = df.T

    # Determine which matplotlib axes to plot on.
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Draw the heatmap.
    sns.heatmap(
        df, ax=ax, vmin=vmin, vmax=vmax, cbar=cbar, cbar_kws=cbar_kws,
        cbar_ax=cbar_ax, square=square, xticklabels=xticklabels,
        yticklabels=yticklabels, **kwargs
    )

    # Update taxa labels, if necessary.
    if pretty_taxa:
        if pname_kws is None:
            pname_kws = {}
        if flip:
            ax.set_yticklabels([common.pname(x, **pname_kws) for x in df.index])
        else:
            ax.set_xticklabels([common.pname(x, **pname_kws) for x in df.columns])

    # Update sample labels, if necessary.
    if label_columns is not None:
        if flip:
            mf = mf.loc[df.columns]
        else:
            mf = mf.loc[df.index]
        mf['sample-id'] = mf.index
        f = lambda r: ' : '.join(r.values.astype(str))
        sample_labels = mf[label_columns].apply(f, axis=1)
        if flip:
            ax.set_xticklabels(sample_labels)
        else:
            ax.set_yticklabels(sample_labels)

    return ax

def clustermap(
    artifact, metadata=None, flip=False, hue1=None, hue_order1=None,
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
    flip : bool, default: False
        If True, flip the x and y axes.
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
    normalize : {None, 'log10', 'clr', 'zscore'}, default: None
        Whether to normalize the the input feature table:

        * None: Do not normalize.
        * 'log10': Apply the log10 transformation adding a psuedocount of 1.
        * 'clr': Apply the centre log ratio (CLR) transformation adding a psuedocount of 1.
        * 'zscore': Apply the Z score transformation.

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

    See Also
    --------
    dokdo.api.clustermap.heatmap

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
        dokdo.clustermap(qza_file,
                         normalize='log10')

    .. image:: images/clustermap-1.png

    We can color the samples by ``body-site`` and use the centered log-ratio
    transformation (CLR) for normalziation:

    .. code:: python3

        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'
        dokdo.clustermap(qza_file,
                         metadata=metadata_file,
                         normalize='clr',
                         hue1='body-site')

    .. image:: images/clustermap-2.png

    We can omit the clustering of samples:

    .. code:: python3

        dokdo.clustermap(qza_file,
                         metadata=metadata_file,
                         normalize='clr',
                         hue1='body-site',
                         row_cluster=False)

    .. image:: images/clustermap-3.png

    We can add an additional grouping variable ``subject``. This may require
    you to use ``bbox_inches='tight'`` option in the :meth:`plt.savefig`
    method when saving the figure as a file (e.g. PNG) because otherwise
    there could be weird issues with the legends (because there are now two
    legends instead of just one). Additionally, note that ``xticklabels``
    and ``yticklabels`` are extra keyword arguments that are passed to
    :meth:`seaborn.clustermap`.

    .. code:: python3

        dokdo.clustermap(qza_file,
                         metadata=metadata_file,
                         normalize='clr',
                         hue1='body-site',
                         hue2='subject',
                         xticklabels=False,
                         yticklabels=False)
        plt.savefig('out.png', bbox_inches='tight')

    .. image:: images/clustermap-4.png

    Finally, we can provide :class:`pandas.DataFrame` as input:

    .. code:: python3

        import pandas as pd
        csv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/table.csv'
        df = pd.read_csv(csv_file, index_col=0)
        dokdo.clustermap(df,
                         metadata=metadata_file,
                         normalize='clr',
                         hue1='body-site')

    .. image:: images/clustermap-5.png
    """
    df = utils.import_feature_table(artifact)

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

    if normalize is not None:
        df = utils.normalize_feature_table(df, normalize)

    # Flip the axes.
    if flip:
        df = df.T
        g = sns.clustermap(df, method=method, metric=metric, figsize=figsize,
                           row_cluster=row_cluster, col_cluster=col_cluster,
                           col_colors=row_colors, **kwargs)
    else:
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
