import tempfile

from . import common

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def taxa_cols(df):
    """Returns metadata columns from DataFrame object."""
    cols = []
    for col in df.columns:
        if 'Unassigned' in col:
            cols.append(col)
        elif '__' in col:
            cols.append(col)
        else:
            continue
    return cols

def _get_mf_cols(df):
    """Returns metadata columns from DataFrame object."""
    cols = []
    for column in df.columns:
        if 'Unassigned' in column:
            continue
        elif '__' in column:
            continue
        else:
            cols.append(column)
    return cols

def _filter_samples(df, mf, exclude_samples, include_samples):
    """Returns DataFrame objects after sample filtering."""
    if exclude_samples and include_samples:
        m = ("Cannot use 'exclude_samples' and "
             "'include_samples' arguments together")
        raise ValueError(m)
    elif exclude_samples:
        for x in exclude_samples:
            for y in exclude_samples[x]:
                i = mf[x] != y
                df = df.loc[i]
                mf = mf.loc[i]
    elif include_samples:
        for x in include_samples:
            i = mf[x].isin(include_samples[x])
            df = df.loc[i]
            mf = mf.loc[i]
    else:
        pass
    return (df, mf)

def _sort_by_mean(df):
    """Returns DataFrame object after sorting taxa by mean relative abundance."""
    a = df.div(df.sum(axis=1), axis=0)
    a = a.loc[:, a.mean().sort_values(ascending=False).index]
    return df[a.columns]

def _get_others_col(df, count, taxa_names, show_others):
    """Returns DataFrame object after selecting taxa."""
    if count != 0 and taxa_names is not None:
        m = "Cannot use 'count' and 'taxa_names' arguments together"
        raise ValueError(m)
    elif count != 0:
        if count < df.shape[1]:
            others = df.iloc[:, count-1:].sum(axis=1)
            df = df.iloc[:, :count-1]
            if show_others:
                df = df.assign(Others=others)
        else:
            pass
    elif taxa_names is not None:
        others = df.drop(columns=taxa_names).sum(axis=1)
        df = df[taxa_names]
        if show_others:
            df = df.assign(Others=others)
    else:
        pass

    return df

def taxa_abundance_bar_plot(
    visualization, metadata=None, level=1, group=None, group_order=None, by=None,
    ax=None, figsize=None, width=0.8, count=0, exclude_samples=None,
    include_samples=None, exclude_taxa=None, sort_by_names=False,
    colors=None, label_columns=None, orders=None, sample_names=None,
    csv_file=None, taxa_names=None, sort_by_mean1=True,
    sort_by_mean2=True, sort_by_mean3=True, show_others=True,
    cmap_name='Accent', legend_short=False, pname_kws=None, legend=True
):
    """
    Create a bar plot showing relative taxa abundance for individual samples.

    The input visualization may already contain sample metadata. To provide
    new sample metadata, and ignore the existing one, use the ``metadata``
    option.

    By default, the method will draw a bar for each sample. To plot the
    average taxa abundance of each sample group, use the ``group`` option.

    .. warning::
        You may get unexpected results when using ``group`` if samples have
        uneven sequencing depth. For example, an outlier sample with an
        extraordinarily large depth could easily mask the contributions from
        the rest of the samples. Therefore, it's strongly recommended to
        rarefy the input feature table when plotting with ``group``.

    +----------------+-----------------------------------------------------+
    | q2-taxa plugin | Example                                             |
    +================+=====================================================+
    | QIIME 2 CLI    | qiime taxa barplot [OPTIONS]                        |
    +----------------+-----------------------------------------------------+
    | QIIME 2 API    | from qiime2.plugins.taxa.visualizers import barplot |
    +----------------+-----------------------------------------------------+

    Parameters
    ----------
    visualization : str, qiime2.Visualization, or pandas.DataFrame
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
        List of metadata columns to be concatenated to form new sample
        labels. Use the string 'sample-id' to indicate the sample ID column.
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
    pname_kws : dict, optional
        Keyword arguments for :meth:`dokdo.api.pname` when ``legend_short``
        is True.
    legend : bool, default: True
        Whether to plot the legend.

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

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()

        qzv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/taxa-bar-plots.qzv'

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            figsize=(10, 7)
        )

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-1.png

    We can change the taxonomic rank from kingdom to genus by setting
    ``level=6``. Note that we are using ``legend=False`` because
    otherwise there will be too many taxa to display on the legend.
    Note also that the colors are recycled in each bar.

    .. code:: python3

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            figsize=(10, 7),
            level=6,
            legend=False
        )

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-2.png

    We can only show the top seven most abundant genera plus 'Others' with
    ``count=8``.

    .. code:: python3

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            figsize=(10, 7),
            level=6,
            count=8,
            legend_short=True
        )

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-3.png

    We can plot the figure and the legend separately.

    .. code:: python3

        fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(10, 7), gridspec_kw={'width_ratios': [9, 1]})

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            ax=ax1,
            level=6,
            count=8,
            legend=False
        )

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            ax=ax2,
            level=6,
            count=8,
            legend_short=True
        )

        handles, labels = ax2.get_legend_handles_labels()

        ax2.clear()
        ax2.legend(handles, labels)
        ax2.axis('off')

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-4.png

    We can use a different color map to display more unique genera (e.g. 20).

    .. code:: python3

        fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(10, 7), gridspec_kw={'width_ratios': [9, 1]})

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            ax=ax1,
            level=6,
            count=20,
            cmap_name='tab20',
            legend=False
        )

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            ax=ax2,
            level=6,
            count=20,
            cmap_name='tab20',
            legend_short=True
        )

        handles, labels = ax2.get_legend_handles_labels()

        ax2.clear()
        ax2.legend(handles, labels)
        ax2.axis('off')

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-5.png

    We can sort the samples by the body-site column in metadata with
    ``by=['body-site']``. To check whether the sorting worked properly,
    we can change the x-axis tick labels to include each sample's
    body-site with ``label_columns``.

    .. code:: python3

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            by=['body-site'],
            label_columns=['body-site', 'sample-id'],
            figsize=(10, 7),
            level=6,
            count=8,
            legend_short=True
        )

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-6.png

    If you want to sort the samples in a certain order instead of ordering
    numerically or alphabetically, use the ``orders`` option.

    .. code:: python3

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            by=['body-site'],
            label_columns=['body-site', 'sample-id'],
            figsize=(10, 7),
            level=6,
            count=8,
            orders={'body-site': ['left palm', 'tongue', 'gut', 'right palm']},
            legend_short=True
        )

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-7.png

    We can only display the 'gut' and 'tongue' samples with
    ``include_samples``.

    .. code:: python3

        fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(9, 7), gridspec_kw={'width_ratios': [9, 1]})

        kwargs = dict(
            include_samples={'body-site': ['gut', 'tongue']},
            by=['body-site'],
            label_columns=['body-site', 'sample-id'],
            level=6,
            count=8
        )

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            ax=ax1,
            legend=False,
            **kwargs
        )

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            ax=ax2,
            legend_short=True,
            **kwargs
        )

        handles, labels = ax2.get_legend_handles_labels()

        ax2.clear()
        ax2.legend(handles, labels)
        ax2.axis('off')

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-8.png

    We can make multiple bar charts grouped by body-site. When making a
    grouped bar chart, it's important to include ``sort_by_mean2=False``
    in order to have the same bar colors for the same taxa across different
    groups.

    .. code:: python3

        fig, axes = plt.subplots(1, 5, figsize=(16, 7))

        groups = ['gut', 'left palm', 'right palm', 'tongue']
        kwargs = dict(level=6, count=8, sort_by_mean2=False, legend=False)

        for i, group in enumerate(groups):
            dokdo.taxa_abundance_bar_plot(
                qzv_file,
                ax=axes[i],
                include_samples={'body-site': [group]},
                **kwargs
            )
            if i != 0:
                axes[i].set_ylabel('')
                axes[i].set_yticks([])
            axes[i].set_title(group)

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            ax=axes[4],
            legend_short=True,
            **kwargs
        )

        handles, labels = axes[4].get_legend_handles_labels()

        axes[4].clear()
        axes[4].legend(handles, labels, loc='center left')
        axes[4].axis('off')

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-9.png

    We can select specific samples with ``sample_names``.

    .. code:: python3

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            figsize=(10, 7),
            level=6,
            count=8,
            sample_names=['L2S382', 'L4S112', 'L1S281'],
            legend_short=True
        )

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-10.png

    We can also pick specific colors for the bars.

    .. code:: python3

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            figsize=(10, 7),
            level=6,
            count=8,
            sample_names=['L2S382', 'L4S112', 'L1S281'],
            colors=['tab:blue', 'tab:orange', 'tab:gray'],
            legend_short=True
        )

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-11.png

    We can create a bar for each sample type.

    .. code:: python3

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            level=6,
            count=8,
            group='body-site',
            figsize=(10, 7),
            legend_short=True
        )

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-12.png

    Of course, we can specify which groups to plot on.

    .. code:: python3

        dokdo.taxa_abundance_bar_plot(
            qzv_file,
            level=6,
            count=8,
            group='body-site',
            group_order=['tongue', 'left palm'],
            figsize=(10, 7),
            legend_short=True
        )

        plt.tight_layout()

    .. image:: images/taxa_abundance_bar_plot-13.png
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
        if pname_kws is None:
            pname_kws = {}
        df.columns = [common.pname(x, **pname_kws) for x in df.columns]

    df.plot.bar(
        stacked=True, legend=legend, ax=ax, width=width, color=c, linewidth=0
    )

    ax.set_xlabel('')
    ax.set_ylabel('Relative abundance (%)')
    if label_columns is not None:
        f = lambda row: ' : '.join(row.values.astype(str))
        xticklabels = mf[label_columns].apply(f, axis=1).tolist()
        ax.set_xticklabels(xticklabels)

    return ax

def taxa_abundance_box_plot(
    visualization, metadata=None, hue=None, hue_order=None,
    level=1, by=None, count=0, exclude_samples=None,
    include_samples=None, exclude_taxa=None, sort_by_names=False,
    sample_names=None, csv_file=None, pseudocount=False,
    taxa_names=None, pretty_taxa=False, show_means=False,
    meanprops=None, show_others=True, sort_by_mean=True,
    add_datapoints=False, jitter=1, alpha=None, size=5, palette=None,
    ax=None, figsize=None,
):
    """
    Create a box plot showing the distribution of relative abundance for
    individual taxa.

    The input visualization may already contain sample metadata. To provide
    new sample metadata, and ignore the existing one, use the ``metadata``
    option.

    By default, the method will draw a box for each taxon. To plot grouped
    box plots, use the ``hue`` option.

    +----------------+-----------------------------------------------------+
    | q2-taxa plugin | Example                                             |
    +================+=====================================================+
    | QIIME 2 CLI    | qiime taxa barplot [OPTIONS]                        |
    +----------------+-----------------------------------------------------+
    | QIIME 2 API    | from qiime2.plugins.taxa.visualizers import barplot |
    +----------------+-----------------------------------------------------+

    Parameters
    ----------
    visualization : str or qiime2.Visualization
        Visualization file or object from the q2-taxa plugin.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    hue : str, optional
        Grouping variable that will produce boxes with different colors.
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    level : int, default: 1
        Taxonomic level at which the features should be collapsed.
    by : list, optional
        Column name(s) to be used for sorting the samples. Using 'sample-id'
        will sort the samples by their name, in addition to other column
        name(s) that may have been provided. If multiple items are provided,
        sorting will occur by the order of the items.
    count : int, default: 0
        Number of top taxa to display. When 0, display all.
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
    sample_names : list, optional
        List of sample IDs to be included.
    csv_file : str, optional
        Path of the .csv file to output the dataframe to.
    pseudocount : bool, default: False
        Add pseudocount to remove zeros.
    taxa_names : list, optional
        List of taxa names to be displayed.
    pretty_taxa : bool, default: False
        If true, only display the smallest taxa rank in the x-axis labels.
    show_means : bool, default: False
        Add means to the boxes.
    meanprops : dict, optional
        The meanprops argument as in matplotlib.pyplot.boxplot.
    show_others : bool, default: True
        Include the 'Others' category.
    sort_by_mean : bool, default: True
        Sort taxa by their mean relative abundance after sample filtration.
    add_datapoints : bool, default: False
        Show data points on top of the boxes.
    jitter : float, default: 1
        Ignored when ``add_datapoints=False``. Amount of jitter (only along
        the categorical axis) to apply.
    alpha : float, optional
        Ignored when ``add_datapoints=False``. Proportional opacity of the
        points.
    size : float, default: 5.0
        Ignored when ``add_datapoints=False``. Radius of the markers, in
        points.
    palette : palette name, list, or dict
        Box colors.
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
    dokdo.api.taxa_abundance_bar_plot
    dokdo.api.addpairs

    Examples
    --------
    Below is a simple example:

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()

        qzv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/taxa-bar-plots.qzv'
        dokdo.taxa_abundance_box_plot(
            qzv_file,
            level=2,
            figsize=(8, 7)
        )
        plt.tight_layout()

    .. image:: images/taxa_abundance_box_plot-1.png

    We can prettify the taxa names:

    .. code:: python3

        dokdo.taxa_abundance_box_plot(
            qzv_file,
            level=2,
            pretty_taxa=True,
            figsize=(8, 7)
        )
        plt.tight_layout()

    .. image:: images/taxa_abundance_box_plot-2.png

    We can only display a selected number of top taxa:

    .. code:: python3

        dokdo.taxa_abundance_box_plot(
            qzv_file,
            level=2,
            count=4,
            pretty_taxa=True,
            figsize=(8, 7)
        )
        plt.tight_layout()

    .. image:: images/taxa_abundance_box_plot-3.png

    We can add data points on top of the boxes:

    .. code:: python3

        dokdo.taxa_abundance_box_plot(
            qzv_file,
            level=2,
            count=4,
            pretty_taxa=True,
            add_datapoints=True,
            figsize=(8, 7)
        )
        plt.tight_layout()

    .. image:: images/taxa_abundance_box_plot-6.png

    We can also specify which taxa to plot:

    .. code:: python3

        dokdo.taxa_abundance_box_plot(
            qzv_file,
            level=2,
            pretty_taxa=True,
            taxa_names=['k__Bacteria;p__Firmicutes', 'k__Bacteria;p__Proteobacteria'],
            figsize=(8, 7)
        )
        plt.tight_layout()

    .. image:: images/taxa_abundance_box_plot-4.png

    In some cases, it may be desirable to plot abundance in log scale. We can achieve this with:

    .. code:: python3

        ax = dokdo.taxa_abundance_box_plot(
            qzv_file,
            level=2,
            hue='body-site',
            count=4,
            pseudocount=True,
            figsize=(8, 7)
        )
        ax.set_ylim([0.05, 200])
        ax.set_yscale('log')
        plt.tight_layout()

    .. image:: images/taxa_abundance_box_plot-5.png
    """
    with tempfile.TemporaryDirectory() as t:
        common.export(visualization, t)
        df = pd.read_csv(f'{t}/level-{level}.csv', index_col=0)

    # If provided, update the metadata.
    if metadata is None:
        pass
    else:
        mf = common.get_mf(metadata)
        cols = _get_mf_cols(df)
        df.drop(columns=cols, inplace=True)
        df = pd.concat([df, mf], axis=1, join='inner')

    df["sample-id"] = df.index

    # If provided, sort the samples for display in the x-axis.
    if by:
        df = df.sort_values(by=by)

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

    df, mf = _filter_samples(df, mf, exclude_samples, include_samples)

    # If provided, only include the specified samples.
    if isinstance(sample_names, list):
        df = df.loc[sample_names]
        mf = mf.loc[sample_names]

    if sort_by_mean:
        df = _sort_by_mean(df)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Add a pseudocount.
    if pseudocount:
        df = df + 1

    # Convert counts to proportions.
    df = df.div(df.sum(axis=1), axis=0)

    df = _get_others_col(df, count, taxa_names, show_others)

    if sort_by_names:
        df = df.reindex(sorted(df.columns), axis=1)

    _taxa_names = df.columns

    df = df * 100

    if hue is not None:
        df2 = pd.concat([df, mf[hue]], axis=1, join='inner')
        df2 = pd.melt(df2, id_vars=[hue])
    else:
        df2 = pd.melt(df)

    if meanprops:
        _meanprops = meanprops
    else:
        _meanprops={'marker':'x',
                    'markerfacecolor':'white',
                    'markeredgecolor':'white',
                    'markersize':'10'}

    d = {}

    if show_means:
        d['showmeans'] = True
        d['meanprops'] = _meanprops

    sns.boxplot(
        x='variable', y='value', hue=hue, hue_order=hue_order, data=df2,
        ax=ax, palette=palette, **d
    )

    if add_datapoints:
        sns.stripplot(
            x='variable', y='value', hue=hue, hue_order=hue_order, data=df2,
            ax=ax, color='black', size=size, dodge=True, jitter=jitter,
            alpha=alpha
        )

    # If provided, output the dataframe as a .csv file.
    if csv_file is not None:
        df3 = pd.concat([df, mf], axis=1, join='inner')
        df3.to_csv(csv_file)

    if pretty_taxa:
        l = [common.pname(x.get_text()) for x in ax.get_xticklabels()]
        ax.set_xticklabels(l)

    ax.set_xlabel('')
    ax.set_ylabel('Relative abundance (%)')

    for ticklabel in ax.get_xticklabels():
        ticklabel.set_rotation(90)

    return ax
