import tempfile
from .common import (_artist, taxa_cols,
    _get_mf_cols, _filter_samples, _sort_by_mean, _get_others_col)
from . import common

import pandas as pd
import matplotlib.pyplot as plt

def taxa_abundance_bar_plot(
    visualization, metadata=None, level=1, group=None, group_order=None, by=None,
    ax=None, figsize=None, width=0.8, count=0, exclude_samples=None,
    include_samples=None, exclude_taxa=None, sort_by_names=False,
    colors=None, label_columns=None, orders=None, sample_names=None,
    csv_file=None, taxa_names=None, sort_by_mean1=True,
    sort_by_mean2=True, sort_by_mean3=True, show_others=True,
    cmap_name='Accent', legend_short=False, artist_kwargs=None
):
    """
    Create a bar plot showing relative taxa abundance for individual samples.

    The input visualization may already contain sample metadata. To provide
    new sample metadata, and ignore the existing one, use the ``metadata``
    option.

    By default, the method will draw a bar for each sample. To plot the
    average taxa abundance of each sample group, use the ``group`` option.

    +----------------+-----------------------------------------------------+
    | q2-taxa plugin | Example                                             |
    +================+=====================================================+
    | QIIME 2 CLI    | qiime taxa barplot [OPTIONS]                        |
    +----------------+-----------------------------------------------------+
    | QIIME 2 API    | from qiime2.plugins.taxa.visualizers import barplot |
    +----------------+-----------------------------------------------------+

    Parameters
    ----------
    visualization : str, qiime2.Visualization, pandas.DataFrame
        Visualization file or object from the q2-taxa plugin. Alternatively,
        a :class:`pandas.DataFrame` object.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    level : int, default: 1
        Taxonomic level at which the features should be collapsed.
    group : str, optional
        Metadata column to be used for grouping the samples.
    group_order : list, optional
        Order to plot the groups in.
    by : list, optional
        Column name(s) to be used for sorting the samples. Using 'sample-id'
        will sort the samples by their name, in addition to other column
        name(s) that may have been provided. If multiple items are provided,
        sorting will occur by the order of the items.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    width : float, default: 0.8
        The width of the bars.
    count : int, default: 0
        The number of taxa to display. When 0, display all.
    exclude_samples : dict, optional
        Filtering logic used for sample exclusion.
        Format: {'col': ['item', ...], ...}.
    include_samples : dict, optional
        Filtering logic used for sample inclusion.
        Format: {'col': ['item', ...], ...}.
    exclude_taxa : list, optional
        The taxa names to be excluded when matched. Case insenstivie.
    sort_by_names : bool, default: False
        If true, sort the columns (i.e. species) to be displayed by name.
    colors : list, optional
        The bar colors.
    label_columns : list, optional
        The column names to be used as the x-axis labels.
    orders : dict, optional
        Dictionary of {column1: [element1, element2, ...], column2:
        [element1, element2...], ...} to indicate the order of items. Used to
        sort the sampels by the user-specified order instead of ordering
        numerically or alphabetically.
    sample_names : list, optional
        List of sample IDs to be included.
    csv_file : str, optional
        Path of the .csv file to output the dataframe to.
    taxa_names : list, optional
        List of taxa names to be displayed.
    sort_by_mean1 : bool, default: True
        Sort taxa by their mean relative abundance before sample filtration.
    sort_by_mean2 : bool, default: True
        Sort taxa by their mean relative abundance after sample filtration by
        'include_samples' or 'exclude_samples'.
    sort_by_mean3 : bool, default: True
        Sort taxa by their mean relative abundance after sample filtration by
        'sample_names'.
    show_others : bool, default: True
        Include the 'Others' category.
    cmap_name : str, default: 'Accent'
        Name of the colormap passed to `matplotlib.cm.get_cmap()`.
    legend_short : bool, default: False
        If true, only display the smallest taxa rank in the legend.
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    dokdo.api.taxa_abundance_box_plot

    Examples
    --------
    Below is a simple example showing taxonomic abundance at the kingdom
    level (i.e. ``level=1``), which is the default taxonomic rank.

    >>> qzv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/taxa-bar-plots.qzv'
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               figsize=(10, 7),
    ...                               artist_kwargs=dict(show_legend=True))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-1.png

    We can change the taxonomic rank from kingdom to genus by setting
    ``level=6``. Note that I removed ``show_legend=True`` because
    otherwise there will be too many taxa to display on the legend.
    Note also that the colors are recycled in each bar.

    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               figsize=(10, 7),
    ...                               level=6)
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-2.png

    We can only show the top seven most abundant genera plus 'Others' with
    ``count=8``.

    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               figsize=(10, 7),
    ...                               level=6,
    ...                               count=8,
    ...                               legend_short=True,
    ...                               artist_kwargs=dict(show_legend=True,
    ...                                                  legend_loc='upper left'))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-3.png

    We can plot the figure and the legend separately.

    >>> fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7), gridspec_kw={'width_ratios': [9, 1]})
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax1,
    ...                               level=6,
    ...                               count=8)
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax2,
    ...                               level=6,
    ...                               count=8,
    ...                               legend_short=True,
    ...                               artist_kwargs=dict(legend_only=True,
    ...                                                  legend_loc='upper left'))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-4.png

    We can use a different color map to display more unique genera (e.g. 20).

    >>> fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(12, 7), gridspec_kw={'width_ratios': [9, 1]})
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax1,
    ...                               level=6,
    ...                               count=20,
    ...                               cmap_name='tab20')
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax2,
    ...                               level=6,
    ...                               count=20,
    ...                               cmap_name='tab20',
    ...                               legend_short=True,
    ...                               artist_kwargs=dict(legend_only=True,
    ...                                                  legend_loc='upper left'))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-5.png

    We can sort the samples by the body-site column in metadata with
    ``by=['body-site']``. To check whether the sorting worked properly,
    we can change the x-axis tick labels to include each sample's
    body-site with ``label_columns``.

    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               by=['body-site'],
    ...                               label_columns=['body-site', 'sample-id'],
    ...                               figsize=(10, 7),
    ...                               level=6,
    ...                               count=8,
    ...                               legend_short=True,
    ...                               artist_kwargs=dict(show_legend=True,
    ...                                                  legend_loc='upper left'))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-6.png

    If you want to sort the samples in a certain order instead of ordering
    numerically or alphabetically, use the ``orders`` option.

    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               by=['body-site'],
    ...                               label_columns=['body-site', 'sample-id'],
    ...                               figsize=(10, 7),
    ...                               level=6,
    ...                               count=8,
    ...                               orders={'body-site': ['left palm', 'tongue', 'gut', 'right palm']},
    ...                               legend_short=True,
    ...                               artist_kwargs=dict(show_legend=True,
    ...                                                  legend_loc='upper left'))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-7.png

    We can only display the 'gut' and 'tongue' samples with
    ``include_samples``.

    >>> fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(9, 7), gridspec_kw={'width_ratios': [9, 1]})
    >>> kwargs = dict(include_samples={'body-site': ['gut', 'tongue']},
    ...                 by=['body-site'],
    ...                 label_columns=['body-site', 'sample-id'],
    ...                 level=6,
    ...                 count=8)
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax1,
    ...                               **kwargs)
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax2,
    ...                               **kwargs,
    ...                               legend_short=True,
    ...                               artist_kwargs=dict(legend_only=True,
    ...                                                  legend_loc='upper left'))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-8.png

    We can make multiple bar charts grouped by body-site. When making a
    grouped bar chart, it's important to include ``sort_by_mean2=False``
    in order to have the same bar colors for the same taxa across different
    groups.

    >>> fig, [ax1, ax2, ax3, ax4, ax5] = plt.subplots(1, 5, figsize=(16, 7), gridspec_kw={'width_ratios': [2, 2, 2, 2, 1]})
    >>> kwargs = dict(level=6, count=8, sort_by_mean2=False)
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax1,
    ...                               include_samples={'body-site': ['gut']},
    ...                               **kwargs,
    ...                               artist_kwargs=dict(title='gut'))
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax2,
    ...                               include_samples={'body-site': ['left palm']},
    ...                               **kwargs,
    ...                               artist_kwargs=dict(title='left palm',
    ...                                                  hide_ylabel=True,
    ...                                                  hide_yticks=True))
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax3,
    ...                               include_samples={'body-site': ['right palm']},
    ...                               **kwargs,
    ...                               artist_kwargs=dict(title='right palm',
    ...                                                  hide_ylabel=True,
    ...                                                  hide_yticks=True))
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax4,
    ...                               include_samples={'body-site': ['tongue']},
    ...                               **kwargs,
    ...                               artist_kwargs=dict(title='tongue',
    ...                                                  hide_ylabel=True,
    ...                                                  hide_yticks=True))
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax5,
    ...                               **kwargs,
    ...                               legend_short=True,
    ...                               artist_kwargs=dict(legend_only=True,
    ...                                                  legend_loc='upper left'))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-9.png

    We can select specific samples with ``sample_names``. We can also
    manually set the x-axis tick labels with ``xticklabels``. Finally, you
    can pick specific colors for the bars.

    >>> fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(10, 5))
    >>> kwargs = dict(level=6, count=3, legend_short=True, sample_names=['L2S382', 'L4S112'])
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax1,
    ...                               **kwargs,
    ...                               artist_kwargs=dict(show_legend=True,
    ...                                                  legend_loc='upper right',
    ...                                                  title="sample_names=['L2S382', 'L4S112']"))
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax2,
    ...                               **kwargs,
    ...                               artist_kwargs=dict(show_legend=True,
    ...                                                  legend_loc='upper right',
    ...                                                  title="xticklabels=['A', 'B']",
    ...                                                  xticklabels=['A', 'B']))
    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               ax=ax3,
    ...                               colors=['tab:blue', 'tab:orange', 'tab:gray'],
    ...                               **kwargs,
    ...                               artist_kwargs=dict(show_legend=True,
    ...                                                  legend_loc='upper right',
    ...                                                  title="colors=['tab:blue', 'tab:orange', 'tab:gray']"))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-10.png

    Finally, we can create a bar for each sample type.

    >>> dokdo.taxa_abundance_bar_plot(qzv_file,
    ...                               level=6,
    ...                               count=8,
    ...                               group='body-site',
    ...                               figsize=(10, 7),
    ...                               legend_short=True,
    ...                               artist_kwargs=dict(show_legend=True))
    >>> plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-11.png
    """
    if isinstance(visualization, pd.DataFrame):
        df = visualization
    else:
        with tempfile.TemporaryDirectory() as t:
            common.export(visualization, t)
            df = pd.read_csv(f'{t}/level-{level}.csv', index_col=0)

    if sort_by_mean1:
        cols = _get_mf_cols(df)
        mf = df[cols]
        df = df.drop(columns=cols)
        df = _sort_by_mean(df)
        df = pd.concat([df, mf], axis=1, join='inner')

    # If provided, update the metadata.
    if metadata is None:
        pass
    else:
        mf = common.get_mf(metadata)
        cols = _get_mf_cols(df)
        df.drop(columns=cols, inplace=True)
        df = pd.concat([df, mf], axis=1, join='inner')

    # If provided, sort the samples by the user-specified order instead of
    # ordering numerically or alphabetically. To do this, we will first add a
    # new temporary column filled with the indicies of the user-provided
    # list. This column will be used for sorting the samples later instead of
    # the original column. After sorting, the new column will be dropped from
    # the dataframe and the original column will replace its place.
    if isinstance(orders, dict):
        for k, v in orders.items():
            u = df[k].unique().tolist()

            if set(u) != set(v):
                message = (f"Target values {u} not matched with user-provided "
                           f"values {v} for metadata column `{k}`")
                raise ValueError(message)

            l = [x for x in range(len(v))]
            d = dict(zip(v, l))
            df.rename(columns={k: f'@{k}'}, inplace=True)
            df[k] = df[f'@{k}'].map(d)

    df["sample-id"] = df.index

    # If provided, sort the samples for display in the x-axis.
    if isinstance(by, list):
        df = df.sort_values(by=by)

    # If sorting was performed by the user-specified order, remove the
    # temporary columns and then bring back the original column.
    if isinstance(orders, dict):
        for k in orders:
            df.drop(columns=[k], inplace=True)
            df.rename(columns={f'@{k}': k}, inplace=True)

    # If provided, exclude the specified taxa.
    if isinstance(exclude_taxa, list):
        dropped = []
        for tax in exclude_taxa:
            for col in df.columns:
                if tax.lower() in col.lower():
                    dropped.append(col)
        dropped = list(set(dropped))
        df = df.drop(columns=dropped)

    # If provided, group the samples by the given metadata column.
    if group is not None:
        df = df.groupby(group)[taxa_cols(df)].agg('sum')

    # Remove the metadata columns.
    cols = _get_mf_cols(df)
    mf = df[cols]
    df = df.drop(columns=cols)

    if group is not None and group_order is not None:
        df = df.loc[group_order]

    df, mf = _filter_samples(df, mf, exclude_samples, include_samples)

    if sort_by_mean2:
        df = _sort_by_mean(df)

    # If provided, only include the specified samples.
    if isinstance(sample_names, list):
        df = df.loc[sample_names]
        mf = mf.loc[sample_names]

        if sort_by_mean3:
            df = _sort_by_mean(df)

    # Convert counts to proportions.
    df = df.div(df.sum(axis=1), axis=0)

    df = _get_others_col(df, count, taxa_names, show_others)

    if sort_by_names:
        df = df.reindex(sorted(df.columns), axis=1)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if isinstance(colors, list):
        c = colors
    else:
        c = plt.cm.get_cmap(cmap_name).colors

    df = df * 100

    # If provided, output the dataframe as a .csv file.
    if csv_file is not None:
        df.to_csv(csv_file)

    if legend_short:
        df.columns = [common.pname(x) for x in df.columns]

    df.plot.bar(stacked=True,
                legend=False,
                ax=ax,
                width=width,
                color=c,
                linewidth=0)

    if label_columns is not None:
        f = lambda row: ' : '.join(row.values.astype(str))
        xticklabels = mf[label_columns].apply(f, axis=1).tolist()
    else:
        xticklabels = None

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': '',
                     'ylabel': 'Relative abundance (%)',
                     'xticklabels': xticklabels,
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
