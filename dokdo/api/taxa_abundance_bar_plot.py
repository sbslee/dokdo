import tempfile
from .common import (_parse_input, _artist,
    _get_mf_cols, _filter_samples, _sort_by_mean, _get_others_col)
import dokdo
import pandas as pd
import matplotlib.pyplot as plt

def taxa_abundance_bar_plot(taxa,
                            metadata=None,
                            level=1,
                            by=None,
                            ax=None,
                            figsize=None,
                            width=0.8,
                            count=0,
                            exclude_samples=None,
                            include_samples=None,
                            exclude_taxa=None,
                            sort_by_names=False,
                            colors=None,
                            label_columns=None,
                            orders=None,
                            sample_names=None,
                            csv_file=None,
                            taxa_names=None,
                            sort_by_mean1=True,
                            sort_by_mean2=True,
                            sort_by_mean3=True,
                            show_others=True,
                            cmap_name='Accent',
                            legend_short=False,
                            artist_kwargs=None):
    """Create a taxa abundance plot.

    Although the input visualization file should contain medatadata already,
    you can replace it with new metadata by using the 'metadata' option.

    Parameters
    ----------
    taxa : str or qiime2.Visualization
        Visualization file or object from the q2-taxa plugin.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    level : int, default: 1
        Taxonomic level at which the features should be collapsed.
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
    taxa_abundance_box_plot

    Notes
    -----
    Example usage of the q2-taxa plugin:
        CLI -> qiime taxa barplot [OPTIONS]
        API -> from qiime2.plugins.taxa.visualizers import barplot

    Examples
    --------
    Plot a heatmap for a numpy array:

    .. plot::
        :context: close-figs

        >>> import numpy as np; np.random.seed(0)
        >>> import seaborn as sns
        >>> uniform_data = np.random.rand(10, 12)
        >>> ax = sns.heatmap(uniform_data)
    """
    with tempfile.TemporaryDirectory() as t:
        _parse_input(taxa, t)
        df = pd.read_csv(f'{t}/level-{level}.csv', index_col=0)

    # If provided, update the metadata.
    if metadata is None:
        pass
    else:
        mf = dokdo.get_mf(metadata)
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

    # Remove the metadata columns.
    cols = _get_mf_cols(df)
    mf = df[cols]
    df = df.drop(columns=cols)

    if sort_by_mean1:
        df = _sort_by_mean(df)

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
        df.columns = [dokdo.pname(x) for x in df.columns]

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
