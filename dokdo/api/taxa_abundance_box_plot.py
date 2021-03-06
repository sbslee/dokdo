import tempfile
from .common import (_parse_input, _artist, _pretty_taxa,
    _get_mf_cols, _filter_samples, _sort_by_mean, _get_others_col)
import dokdo
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def taxa_abundance_box_plot(taxa,
                            metadata=None,
                            hue=None,
                            hue_order=None,
                            add_datapoints=False,
                            level=1,
                            by=None,
                            ax=None,
                            figsize=None,
                            count=0,
                            exclude_samples=None,
                            include_samples=None,
                            exclude_taxa=None,
                            sort_by_names=False,
                            sample_names=None,
                            csv_file=None,
                            size=5,
                            pseudocount=False,
                            taxa_names=None,
                            brief_xlabels=False,
                            show_means=False,
                            meanprops=None,
                            show_others=True,
                            sort_by_mean=True,
                            jitter=1,
                            alpha=None,
                            artist_kwargs=None):
    """
    This method creates a taxa abundance box plot.

    Parameters
    ----------
    taxa : str or qiime2.Visualization
        Visualization file or object from the q2-taxa plugin.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    hue : str, optional
        Grouping variable that will produce boxes with different colors.
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    add_datapoints : bool, default: False
        Show datapoints on top of the boxes.
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
    sample_names : list, optional
        List of sample IDs to be included.
    csv_file : str, optional
        Path of the .csv file to output the dataframe to.
    size : float, default: 5.0
        Radius of the markers, in points.
    pseudocount : bool, default: False
        Add pseudocount to remove zeros.
    taxa_names : list, optional
        List of taxa names to be displayed.
    brief_xlabels : bool, default: False
        If true, only display the smallest taxa rank in the x-axis labels.
    show_means : bool, default: False
        Add means to the boxes.
    meanprops : dict, optional
        The meanprops argument as in matplotlib.pyplot.boxplot.
    show_others : bool, default: True
        Include the 'Others' category.
    sort_by_mean : bool, default: True
        Sort taxa by their mean relative abundance after sample filtration.
    jitter : float, default: 1
        Amount of jitter (only along the categorical axis) to apply.
    alpha : float, optional
        Proportional opacity of the points.
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    taxa_abundance_bar_plot
    addpairs

    Notes
    -----
    Example usage of the q2-taxa plugin:
        CLI -> qiime taxa barplot [OPTIONS]
        API -> from qiime2.plugins.taxa.visualizers import barplot
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

    sns.boxplot(x='variable',
                y='value',
                hue=hue,
                hue_order=hue_order,
                data=df2,
                ax=ax,
                **d)

    if add_datapoints:
        remove_duplicates = True
        # Alternative method: sns.swarmplot()
        sns.stripplot(x='variable',
                      y='value',
                      hue=hue,
                      hue_order=hue_order,
                      data=df2,
                      ax=ax,
                      color='black',
                      size=size,
                      dodge=True,
                      jitter=jitter,
                      alpha=alpha)
    else:
        remove_duplicates = False

    # If provided, output the dataframe as a .csv file.
    if csv_file is not None:
        df3 = pd.concat([df, mf], axis=1, join='inner')
        df3.to_csv(csv_file)

    if brief_xlabels:
        xticklabels = [_pretty_taxa(x) for x in ax.get_xticklabels()]
    else:
        xticklabels = None

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xrot': 45,
                     'xha': 'right',
                     'xlabel': '',
                     'ylabel': 'Relative abundance (%)',
                     'xticklabels': xticklabels,
                     'remove_duplicates': remove_duplicates,
                     **artist_kwargs}

    if hue is not None:
        artist_kwargs['legend_title'] = hue

    ax = _artist(ax, **artist_kwargs)

    return ax
