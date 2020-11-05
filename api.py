# Import standard libraries.
import math
from tempfile import TemporaryDirectory
import warnings

# Import external libraries.
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import skbio as sb

# Import QIIME 2 libraries
import qiime2
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2 import Visualization
from qiime2.plugins import feature_table
from qiime2.plugins import diversity_lib
from qiime2.plugins import diversity










# -- Private methods ---------------------------------------------------------

def _get_mf_cols(df):
    "Returns metadata columns from DataFrame object."
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
    "Returns DataFrame objects after sample filtering."
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
    "Returns DataFrame object after sorting columns by their mean."
    return df.loc[:, df.mean().sort_values(ascending=False).index]










def _pretty_taxa(s):
    "Returns pretty taxa name."
    if isinstance(s, matplotlib.text.Text):
        s = s.get_text()
    ranks = list(reversed(s.split(';')))

    for i, rank in enumerate(ranks):
        if rank in ['Others', 'Unassigned']:
            return rank

        if rank == '__':
            continue

        if rank.split('__')[1] is '':
            continue

        if 'uncultured' in rank:
            continue

        # The species name can be sometimes tricky to parse because it could 
        # be full (e.g. Helicobacter pylori) or partial (e.g. pylori). In the 
        # latter case, I will borrow the genus name (e.g. Helicobacter) to 
        # form the full species name.
        if 's__' in rank:
            rank = rank.split('__')[1]

            if len(rank.split('_')) == 1:
                genus = ranks[i+1].split('__')[1].split('_')[0]
                species = rank.split('_')[0]
                rank = f'{genus} {species}'
            else:
                rank = rank.replace('_', ' ')

        if '__' in rank:
            rank = rank.split('__')[1]

        return rank










def _legend_handler(ax,
                    show_legend,
                    legend_loc,
                    legend_ncol=1,
                    legend_labels=None,
                    legend_short=False,
                    remove_duplicates=False,
                    legend_only=False):
    "Handles the legend of a figure."
    h, l = ax.get_legend_handles_labels()

    if legend_short:
        l = [_pretty_taxa(x) for x in l]

    if legend_labels:
        a = len(legend_labels)
        b = len(l)
        if a != b:
            m = f"Expected {b} legend labels, received {a}"
            raise ValueError(m)
        l = legend_labels

    if remove_duplicates:
        if h:
            n = int(len(h) / 2)
            h, l = h[:n], l[:n]

    if legend_only:
        ax.clear()
        ax.legend(h, l, loc=legend_loc, ncol=legend_ncol)
        ax.axis('off')
    elif show_legend:
        if h:
            ax.legend(h, l, loc=legend_loc, ncol=legend_ncol)
        else:
            warnings.warn("No handles with labels found to put in legend.")
    else:
        if ax.get_legend():
            ax.get_legend().remove()
        else:
            pass

    return ax










def _artist(ax,
            title=None,
            hide_xlabel=False,
            hide_ylabel=False,
            hide_xticks=False,
            hide_yticks=False,
            hide_xticklabels=False,
            hide_yticklabels=False,
            xmin=None,
            xmax=None,
            ymin=None,
            ymax=None,
            xlog=False,
            ylog=False,
            **kwargs):
    """
    This method sets various properties of the given figure.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object to draw the plot onto.
    title : str, optional
        Sets the figure title.
    hide_xlabel : bool, default: False
        Hides the x-axis label.
    hide_ylabel : bool, default: False
        Hides the y-axis label.
    hide_xticks : bool, default: False
        Hides ticks and tick labels for the x-axis.
    hide_yticks : bool, default: False
        Hides ticks and tick labels for the y-axis.
    hide_xticklabels : bool, default: False
        Hides tick labels for the x-axis.
    hide_yticklabels : bool, default: False
        Hides tick labels for the y-axis.
    xmin : float, optional
        Minimum value for the x-axis.
    xmax : float, optional
        Maximum value for the x-axis.
    ymin : float, optional
        Minimum value for the y-axis.
    ymax : float, optional
        Maximum value for the x-axis.
    xlog : bool, default: False
        Draw the x-axis in log scale.
    ylog : bool, default: False
        Draw the y-axis in log scale.

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    if isinstance(title, str):
        ax.set_title(title)

    if hide_xlabel:
        ax.set_xlabel('')

    if hide_ylabel:
        ax.set_ylabel('')

    if hide_xticks:
        ax.set_xticks([])

    if hide_yticks:
        ax.set_yticks([])

    if hide_xticklabels:
        ax.set_xticklabels([])

    if hide_yticklabels:
        ax.set_yticklabels([])

    ax.set_xlim(left=xmin, right=xmax)
    ax.set_ylim(bottom=ymin, top=ymax)

    if xlog:
        ax.set_xscale('log')

    if ylog:
        ax.set_yscale('log')

    return ax










# -- General methods ---------------------------------------------------------

def get_mf(metadata):
    """
    This method automatically detects the type of input metadata and converts 
    it to DataFrame object.

    Parameters
    ----------
    metadata : str or qiime2.metadata.metadata.Metadata
        Metadata file or object.

    Returns
    -------
    pandas.DataFrame
        DataFrame object containing metadata.
    """
    if isinstance(metadata, str):
        mf = Metadata.load(metadata).to_dataframe()
    elif isinstance(metadata, qiime2.metadata.metadata.Metadata):
        mf = metadata.to_dataframe()
    else:
        raise TypeError(f"Incorrect metadata type: {type(metadata)}")
    return mf










def ordinate(table,
             metadata=None,
             where=None,
             metric='jaccard',
             phylogeny=None):
    """
    This method wraps multiple QIIME 2 methods to perform ordination and 
    returns Artifact object containing PCoA results.

    Under the hood, this method performs sample filtration, rarefaction, 
    distance matrix computation, and PCoA analysis.

    Parameters
    ----------
    table : str
        Table file.
    metadata : str, optional
        Metadata file.
    where : str, optional
        SQLite WHERE clause specifying sample metadata criteria.
    metric : str, default: 'jaccard'
        Metric used for distance matrix computation ('jaccard',
        'bray_curtis', 'unweighted_unifrac', or 'weighted_unifrac').
    phylogeny : str, optional
        Rooted tree file. Required if using 'unweighted_unifrac', or 
        'weighted_unifrac' as metric.

    Returns
    -------
    qiime2.sdk.result.Artifact
        Artifact containing PCoA results from 'diversity.methods.pcoa'.
    """
    if where:
        if metadata is None:
            m = "To use 'where' argument, you must provide metadata"
            raise ValueError(m)

        filter_result = feature_table.methods.filter_samples(
            table=Artifact.load(table),
            metadata=Metadata.load(metadata),
            where=where,
        )
        _table = filter_result.filtered_table
    else:
        _table = Artifact.load(table)

    min_depth = int(_table.view(pd.DataFrame).sum(axis=1).min())

    rarefy_result = feature_table.methods.rarefy(table=_table,
                                                 sampling_depth=min_depth)
    
    rarefied_table = rarefy_result.rarefied_table

    if metric == 'jaccard':
        distance_matrix_result = diversity_lib.methods.jaccard(table=rarefied_table)
    elif metric == 'bray_curtis':
        distance_matrix_result = diversity_lib.methods.bray_curtis(table=rarefied_table)
    elif metric == 'unweighted_unifrac':
        distance_matrix_result = diversity_lib.methods.unweighted_unifrac(table=rarefied_table, phylogeny=Artifact.load(phylogeny))
    elif metric == 'weighted_unifrac':
        distance_matrix_result = diversity_lib.methods.weighted_unifrac(table=rarefied_table, phylogeny=Artifact.load(phylogeny))
    else:
        raise ValueError(f"Incorrect metric detected: {metric}")
    
    distance_matrix = distance_matrix_result.distance_matrix
    
    pcoa_result = diversity.methods.pcoa(distance_matrix=distance_matrix)

    return pcoa_result.pcoa










# -- Main plotting methods ---------------------------------------------------

def read_quality_plot(demux,
                      strand='forward',
                      ax=None,
                      figsize=None,
                      **kwargs):
    """
    This method creates a read quality plot.

    Parameters
    ----------
    demux : str
        Path to the visualization file from the 'qiime demux summarize' 
        command.
    strand : str, default: 'forward'
        Read strand to be displayed (either 'forward' or 'reverse').
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    t = TemporaryDirectory()
    Visualization.load(demux).export_data(t.name)

    l = ['forward', 'reverse']
    if strand not in l:
        raise ValueError(f"Strand should be one of the following: {l}")

    df = pd.read_table(f'{t.name}/{strand}-seven-number-summaries.tsv',
                       index_col=0, skiprows=[1])
    df = pd.melt(df.reset_index(), id_vars=['index'])
    df['variable'] = df['variable'].astype('int64')
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.boxplot(x='variable', y='value', data=df, ax=ax, fliersize=0,
                boxprops=dict(color='white', edgecolor='black'),
                medianprops=dict(color='red'),
                whiskerprops=dict(linestyle=':'))

    ax.set_ylim([0, 45])
    ax.set_xlabel('Sequence base')
    ax.set_ylabel('Quality score')
    a = np.arange(df['variable'].min(), df['variable'].max(), 20)
    ax.set_xticks(a)
    ax.set_xticklabels(a)

    ax = _artist(ax, **kwargs)

    return ax










def denoising_stats_plot(stats,
                         metadata,
                         where,
                         ax=None,
                         figsize=None,
                         pseudocount=False,
                         order=None,
                         hide_nsizes=False,
                         **kwargs):
    """
    This method creates a grouped box plot using denoising statistics from 
    DADA2 (i.e. the 'qiime dada2 denoise-paired' command).

    Parameters
    ----------
    stats : str
        Path to the denoising-stats.qza file.
    metadata : str or qiime2.metadata.metadata.Metadata
        Metadata file or object.
    where : str
        Column name of the sample metadata.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    pseudocount : bool, default: False
        Add pseudocount to remove zeros.
    order : list, optional
        Order to plot the categorical levels in.
    hide_nsizes : bool, default: False
        Hide sample size from x-axis labels.

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    t = TemporaryDirectory()
    Artifact.load(stats).export_data(t.name)
    df1 = pd.read_table(f'{t.name}/stats.tsv', skiprows=[1], index_col=0)

    mf = get_mf(metadata)

    df2 = pd.concat([df1, mf], axis=1, join='inner')

    a = ['input', 'filtered', 'denoised', 'merged', 'non-chimeric', where]
    df3 = pd.melt(df2[a], id_vars=[where])

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if pseudocount:
        df3['value'] = df3['value'] + 1

    sns.boxplot(x=where,
                y='value',
                data=df3,
                hue='variable',
                ax=ax,
                order=order)

    if hide_nsizes is False:
        nsizes = df2[where].value_counts().to_dict()
        xtexts = [x.get_text() for x in ax.get_xticklabels()]
        xtexts = [f'{x} ({nsizes[x]})' for x in xtexts]
        ax.set_xticklabels(xtexts)

    ax = _artist(ax, **kwargs)

    return ax










def alpha_rarefaction_plot(rarefaction,
                           hue='sample-id',
                           metric='shannon',
                           ax=None,
                           figsize=None,
                           show_legend=False,
                           legend_loc='best',
                           legend_ncol=1,
                           legend_only=False,
                           hue_order=None,
                           **kwargs):
    """
    This method creates an alpha rarefaction plot.

    Parameters
    ----------
    rarefaction : str or qiime2.sdk.result.Visualization
        Visualization file or object from alpha rarefaction.
    hue : str, default: 'sample-id'
        Grouping variable that will produce lines with different colors.
    metric : str, default: 'shannon'
        Diversity metric ('shannon', 'observed_features', or 'faith_pd').
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    show_legend : bool, default: False
        Show the legend.
    legend_loc : str, default: 'best'
        Legend location specified as in matplotlib.pyplot.legend.
    legend_ncol : int, default: 1
        The number of columns that the legend has.
    legend_only : bool, default: False
        Plot the legend only.
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    kwargs : dict, optional
        Other keyword arguments passed down to matplotlib.axes.Axes.scatter.

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    t = TemporaryDirectory()

    if isinstance(rarefaction, qiime2.sdk.result.Visualization):
        fn = f'{t.name}/alpha-rarefaction.qzv'
        rarefaction.save(fn)
        rarefaction = fn

    Visualization.load(rarefaction).export_data(t.name)

    l = ['observed_features', 'faith_pd', 'shannon']
    if metric not in l:
        raise ValueError(f"Metric should be one of the following: {l}")

    df = pd.read_csv(f'{t.name}/{metric}.csv', index_col=0)

    metadata_columns = [x for x in df.columns if 'iter' not in x]

    df = pd.melt(df.reset_index(), id_vars=['sample-id'] + metadata_columns)

    df['variable'] = df['variable'].str.split('_').str[0].str.replace(
                         'depth-', '').astype(int)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.lineplot(x='variable',
                 y='value',
                 data=df,
                 hue=hue,
                 ax=ax,
                 err_style='bars',
                 sort=False,
                 hue_order=hue_order)

    ax.set_xlabel('Sequencing depth')
    ax.set_ylabel(metric)

    ax = _legend_handler(ax, show_legend, legend_loc, legend_ncol=legend_ncol, legend_only=legend_only)

    ax = _artist(ax, **kwargs)

    return ax










def alpha_diversity_plot(significance,
                         where,
                         ax=None,
                         figsize=None,
                         add_swarmplot=False,
                         order=None,
                         ylabel=None,
                         **kwargs):
    """
    This method creates an alpha diversity plot.

    Parameters
    ----------
    significance : str
        Path to the visualization file from the 'qiime diversity 
        alpha-group-significance' command.
    where : str
        Column name to be used for the x-axis.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    add_swarmplot : bool, default: False
        Add a swarm plot on top of the box plot.
    order : list, optional
        Order to plot the categorical levels in.
    ylabel : str, optional
        Y-axis label.

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    t = TemporaryDirectory()
    Visualization.load(significance).export_data(t.name)
    df = Metadata.load(f'{t.name}/metadata.tsv').to_dataframe()
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    metric = df.columns[-1]

    boxprops = dict(color='white', edgecolor='black')

    d = {'x': where, 'y': metric, 'ax': ax, 'order': order, 'data': df}

    sns.boxplot(boxprops=boxprops, **d)

    if add_swarmplot:
        sns.swarmplot(**d)

    if ylabel:
        ax.set_ylabel(ylabel)

    ax = _artist(ax, **kwargs)

    return ax










def beta_2d_plot(ordination,
                 metadata=None,
                 hue=None,
                 size=None,
                 style=None,
                 s=80,
                 alpha=None,
                 ax=None,
                 figsize=None,
                 show_legend=False,
                 legend_loc='best',
                 hue_order=None,
                 style_order=None,
                 legend_ncol=1,
                 legend_type='brief',
                 **kwargs):
    """
    This method creates a 2D beta diversity plot.

    Parameters
    ----------
    ordination : str or qiime2.sdk.result.Artifact
        Artifact file or object from ordination.
    metadata : str or qiime2.metadata.metadata.Metadata, optional
        Metadata file or object.
    hue : str, optional
        Grouping variable that will produce points with different colors.
    size : str, optional
        Grouping variable that will produce points with different sizes.
    style : str, optional
        Grouping variable that will produce points with different markers.
    s : int, default: 80
        Marker size.
    alpha : float, optional
        Proportional opacity of the points.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    show_legend : bool, default: False
        Show the legend.
    legend_loc : str, default: 'best'
        Legend location specified as in matplotlib.pyplot.legend.
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    style_order : list, optional
        Specify the order of categorical levels of the 'style' semantic.
    legend_ncol : int, default: 1
        The number of columns that the legend has.
    legend_type : str, default: 'brief'
        Legend type as in seaborn.scatterplot ('brief' or 'full').
    kwargs : dict, optional
        Other keyword arguments passed down to matplotlib.axes.Axes.scatter.

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    t = TemporaryDirectory()

    if isinstance(ordination, qiime2.sdk.result.Artifact):
        fn = f'{t.name}/ordination.qza'
        ordination.save(fn)
        ordination = fn

    Artifact.load(ordination).export_data(t.name)

    df1 = pd.read_table(f'{t.name}/ordination.txt', header=None, index_col=0,
                        skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                        skipfooter=4, engine='python', usecols=[0, 1, 2])
    df1.columns = ['A1', 'A2']

    if metadata is None:
        df2 = df1
    else:
        mf = get_mf(metadata)
        df2 = pd.concat([df1, mf], axis=1, join='inner')

    with open(f'{t.name}/ordination.txt') as f:
        v = [round(float(x) * 100, 2) for x in f.readlines()[4].split('\t')]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.scatterplot(data=df2,
                    x='A1',
                    y='A2',
                    hue=hue,
                    hue_order=hue_order,
                    style=style,
                    style_order=style_order,
                    size=size,
                    ax=ax,
                    s=s,
                    alpha=alpha,
                    legend=legend_type)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(f'Axis 1 ({v[0]} %)')
    ax.set_ylabel(f'Axis 2 ({v[1]} %)')

    ax = _legend_handler(ax, show_legend, legend_loc, legend_ncol=legend_ncol)

    ax = _artist(ax, **kwargs)

    return ax










def beta_3d_plot(ordination,
                 metadata,
                 hue=None,
                 azim=-60,
                 elev=30,
                 s=80, 
                 ax=None,
                 figsize=None,
                 show_legend=False,
                 legend_loc='best',
                 hue_order=None,
                 legend_ncol=1,
                 **kwargs):
    """
    This method creates a 3D beta diversity plot.

    Parameters
    ----------
    ordination : str
        Path to the artifact file from ordination (e.g. 
        bray_curtis_pcoa_results.qza).
    metadata : str or qiime2.metadata.metadata.Metadata
        Metadata file or object.
    hue : str, optional
        Grouping variable that will produce points with different colors.
    azim : int, default: -60
        Elevation viewing angle.
    elev : int, default: 30
        Azimuthal viewing angle.
    s : int, default: 80
        Marker size.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    show_legend : bool, default: False
        Show the legend.
    legend_loc : str, default: 'best'
        Legend location specified as in matplotlib.pyplot.legend.
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    legend_ncol : int, default: 1
        The number of columns that the legend has.
    kwargs : dict, optional
        Other keyword arguments passed down to matplotlib.axes.Axes.scatter.

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """

    t = TemporaryDirectory()
    Artifact.load(ordination).export_data(t.name)

    df = pd.read_table(f'{t.name}/ordination.txt', header=None, index_col=0,
                        skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                        skipfooter=4, engine='python')
    df = df.sort_index()

    mf = get_mf(metadata)
    mf = mf.sort_index()
    mf = mf.assign(**{'sample-id': mf.index})

    with open(f'{t.name}/ordination.txt') as f:
        v = [round(float(x) * 100, 2) for x in f.readlines()[4].split('\t')]

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection='3d')

    ax.view_init(azim=azim, elev=elev)
    ax.set_xlabel(f'Axis 1 ({v[0]} %)')
    ax.set_ylabel(f'Axis 2 ({v[1]} %)')
    ax.set_zlabel(f'Axis 3 ({v[2]} %)')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    d = {'s': s}

    if hue is None:
        ax.scatter(df.iloc[:, 0],
                   df.iloc[:, 1],
                   df.iloc[:, 2],
                   **d)
    else:
        if hue_order is None:
            levels = sorted(mf[hue].unique())
        else:
            levels = hue_order

        for c in levels:
            i = mf[hue] == c
            df2 = df.loc[i]
            ax.scatter(df2.iloc[:, 0],
                       df2.iloc[:, 1],
                       df2.iloc[:, 2],
                       label=c,
                       **d)

    ax = _legend_handler(ax, show_legend, legend_loc, legend_ncol=legend_ncol)

    ax = _artist(ax, **kwargs)

    return ax










def distance_matrix_plot(distance_matrix,
                         bins=100,
                         pairs=None,
                         ax=None,
                         figsize=None,
                         **kwargs):
    """
    This method creates a histogram from a distance matrix.

    Parameters
    ----------
    distance_matrix : str or qiime2.sdk.result.Artifact
         Artifact file or object from distance matrix computation.
    bins : int, optional
        Number of bins to be displayed.
    pairs : list, optional
        List of sample pairs to be shown in red vertical lines.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    t = TemporaryDirectory()

    if isinstance(distance_matrix, qiime2.sdk.result.Artifact):
        fn = f'{t.name}/distance-matrix.qza'
        distance_matrix.save(fn)
        distance_matrix = fn

    Artifact.load(distance_matrix).export_data(t.name)
    df = pd.read_table(f'{t.name}/distance-matrix.tsv', index_col=0)
    dist = sb.stats.distance.DistanceMatrix(df, ids=df.columns)
    cdist = dist.condensed_form()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    ax.hist(cdist, bins=bins)

    # https://stackoverflow.com/a/36867493/7481899
    def square_to_condensed(i, j, n):
        assert i != j, "no diagonal elements in condensed matrix"
        if i < j:
            i, j = j, i
        return n*j - j*(j+1)//2 + i - 1 - j

    if pairs:
        idx = []

        for pair in pairs:
            i = square_to_condensed(dist.index(pair[0]), dist.index(pair[1]), len(dist.ids))
            idx.append(cdist[i])

        for i in idx:
            ax.axvline(x=i, c='red')

    ax = _artist(ax, **kwargs)

    return ax










def taxa_abundance_bar_plot(taxa,
                            metadata=None,
                            level=1,
                            by=[],
                            ax=None,
                            figsize=None,
                            width=0.8,
                            count=0,
                            exclude_samples=None,
                            include_samples=None,
                            exclude_taxa=[],
                            show_legend=False,
                            legend_short=False,
                            legend_loc='best',
                            legend_labels=None,
                            legend_only=False,
                            sort_by_names=False,
                            colors=[],
                            label_columns=[],
                            orders={},
                            sample_names=[],
                            csv_file=None,
                            xlabels=None,
                            taxa_names=None,
                            sort_by_mean1=True,
                            sort_by_mean2=True,
                            **kwargs):
    """
    This method creates a taxa abundance plot.

    Although the input visualization file should contain medatadata already, 
    you can replace it with new metadata by using the 'metadata' option.

    Parameters
    ----------
    taxa : str
        Path to the visualization file from the 'qiime taxa barplot'.
    metadata : str or qiime2.metadata.metadata.Metadata
        Metadata file or object.
    level : int
        Taxonomic level at which the features should be collapsed.
    by : list of str
        Column name(s) to be used for sorting the samples. Using 'index' will 
        sort the samples by their name, in addition to other column name(s) 
        that may have been provided. If multiple items are provided, sorting 
        will occur by the order of the items.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    width : float
        The width of the bars.
    count : int, default: 0
        The number of taxa to display. When 0, display all.
    exclude_samples : dict, optional
        Filtering logic used for sample exclusion.
        Format: {'col': ['item', ...], ...}.
    include_samples : dict, optional
        Filtering logic used for sample inclusion.
        Format: {'col': ['item', ...], ...}.
    exclude_taxa : list
        The taxa names to be excluded when matched. Case insenstivie.
    show_legend : bool, default: False
        Show the legend.
    legend_short : bool
        If true, only display the smallest taxa rank in the legend.
    legend_loc : str, default: 'best'
        Legend location specified as in matplotlib.pyplot.legend.
    legend_labels : list, optional
        Legend texts.
    legend_only : bool, default: False
        Plot the legend only.
    sort_by_names : bool
        If true, sort the columns (i.e. species) to be displayed by name.
    colors : list
        The bar colors.
    label_columns : list
        The column names to be used as the x-axis labels.
    orders : dict
        Dictionary of {column1: [element1, element2, ...], column2: 
        [element1, element2...], ...} to indicate the order of items. Used to 
        sort the sampels by the user-specified order instead of ordering 
        numerically or alphabetically.
    sample_names : list
        List of sample IDs to be included.
    csv_file : str
        Path of the .csv file to output the dataframe to.
    xlabels : list, optional
        List of the x-axis labels.
    taxa_names : list, optional
        List of taxa names to be displayed.
    sort_by_mean1 : bool, default: True
        Sort taxa by their mean abundance before sample filtration.
    sort_by_mean2 : bool, default: True
        Sort taxa by their mean abundance after sample filtration.
    kwargs : dict, optional
        Additional keyword arguments are documented in DataFrame.plot.

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    t = TemporaryDirectory()
    Visualization.load(taxa).export_data(t.name)
    df = pd.read_csv(f'{t.name}/level-{level}.csv', index_col=0)

    # If provided, update the metadata.
    if metadata is None:
        pass
    else:
        mf = get_mf(metadata)
        cols = _get_mf_cols(df)
        df.drop(columns=cols, inplace=True)
        df = pd.concat([df, mf], axis=1, join='inner')

    # If provided, sort the samples by the user-specified order instead of 
    # ordering numerically or alphabetically. To do this, we will first add a 
    # new temporary column filled with the indicies of the user-provided 
    # list. This column will be used for sorting the samples later instead of 
    # the original column. After sorting, the new column will be dropped from 
    # the dataframe and the original column will replace its place.
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

    # If provided, sort the samples for display in the x-axis.
    if by:
        df = df.sort_values(by=by)

    # If sorting was performed by the user-specified order, remove the 
    # temporary columns and then bring back the original column.
    for k in orders:
        df.drop(columns=[k], inplace=True)
        df.rename(columns={f'@{k}': k}, inplace=True)

    # If provided, exclude the specified taxa.
    if exclude_taxa:
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
    mf = mf.assign(**{'sample-id': mf.index})
    df = df.drop(columns=cols)

    if sort_by_mean1:
        a = df.div(df.sum(axis=1), axis=0)
        a = _sort_by_mean(a)
        df = df[a.columns]

    df, mf = _filter_samples(df, mf, exclude_samples, include_samples)

    # If provided, only include the specified samples.
    if sample_names:
        df = df.loc[sample_names]
        mf = mf.loc[sample_names]

    # Convert counts to proportions.
    df = df.div(df.sum(axis=1), axis=0)

    if sort_by_mean2:
        df = _sort_by_mean(df)

    # If provided, collapse species to the Others column.
    if count is not 0 and taxa_names is not None:
        m = "Cannot use 'count' and 'taxa_names' arguments together"
        raise ValueError(m)
    elif count is not 0:
        others = df.iloc[:, count-1:].sum(axis=1)
        df = df.iloc[:, :count-1]
        df['Others'] = others
    elif taxa_names is not None:
        others = df.drop(columns=taxa_names).sum(axis=1)
        df = df[taxa_names]
        df['Others'] = others
    else:
        pass

    if sort_by_names:
        df = df.reindex(sorted(df.columns), axis=1)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if colors:
        c = colors
    else:
        c = plt.cm.get_cmap('Accent').colors

    df = df * 100

    df.plot.bar(stacked=True,
                legend=False,
                ax=ax,
                width=width,
                color=c,
                linewidth=0)

    ax.set_xlabel('')
    ax.set_ylabel('Relative abundance (%)')


    if label_columns is not None:
        f = lambda row: ' : '.join(row.values.astype(str))
        smart_xlabels = mf[label_columns].apply(f, axis=1)
    else:
        smart_xlabels = None



    # If provided, output the dataframe as a .csv file.
    if csv_file is not None:
        df.to_csv(csv_file)

    # Manage the x-axis labels.
    if xlabels is not None:
        xtexts = [x.get_text() for x in ax.get_xticklabels()]
        if len(xtexts) != len(xlabels):
            m = f"Expected {len(xtexts)} items, but found {len(xlabels)}"
            raise ValueError(m)
        ax.set_xticklabels(xlabels)

    ax = _legend_handler(ax,
                         show_legend,
                         legend_loc,
                         legend_labels=legend_labels,
                         legend_short=legend_short,
                         legend_only=legend_only)

    ax = _artist(ax, **kwargs)

    return ax










def taxa_abundance_box_plot(taxa,
                            hue=None,
                            hue_order=None,
                            add_datapoints=False,
                            level=1,
                            by=[],
                            ax=None,
                            figsize=None,
                            count=0,
                            exclude_samples=None,
                            include_samples=None,
                            exclude_taxa=[],
                            sort_by_names=False,
                            sample_names=[],
                            csv_file=None,
                            size=5,
                            xlabels=None,
                            pseudocount=False,
                            taxa_names=None,
                            brief_xlabels=False,
                            show_legend=False,
                            legend_loc='best',
                            show_means=False,
                            meanprops=None,
                            **kwargs):
    """
    This method creates a taxa abundance box plot.

    Parameters
    ----------
    taxa : str
        Path to the visualization file from the 'qiime taxa barplot'.
    hue : str, optional
        Grouping variable that will produce boxes with different colors.
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    add_datapoints : bool, default: False
        Show datapoints on top of the boxes.
    level : int
        Taxonomic level at which the features should be collapsed.
    by : list of str
        Column name(s) to be used for sorting the samples. Using 'index' will 
        sort the samples by their name, in addition to other column name(s) 
        that may have been provided. If multiple items are provided, sorting 
        will occur by the order of the items.
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
    exclude_taxa : list
        The taxa names to be excluded when matched. Case insenstivie.
    sort_by_names : bool
        If true, sort the columns (i.e. species) to be displayed by name.
    sample_names : list
        List of sample IDs to be included.
    csv_file : str
        Path of the .csv file to output the dataframe to.
    size : float, default: 5.0
        Radius of the markers, in points.
    xlabels : list, optional
        List of the x-axis labels.
    pseudocount : bool, default: False
        Add pseudocount to remove zeros.
    taxa_names : list, optional
        List of taxa names to be displayed.
    brief_xlabels : bool, default: False
        If true, only display the smallest taxa rank in the x-axis labels.
    show_legend : bool, default: False
        Show the legend.
    legend_loc : str, default: 'best'
        Legend location specified as in matplotlib.pyplot.legend.
    show_means : bool, default: False
        Add means to the boxes.
    meanprops : dict, optional
        The meanprops argument as in matplotlib.pyplot.boxplot.

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    t = TemporaryDirectory()
    Visualization.load(taxa).export_data(t.name)
    df = pd.read_csv(f'{t.name}/level-{level}.csv', index_col=0)

    # If provided, sort the samples for display in the x-axis.
    if by:
        df = df.sort_values(by=by)

    # If provided, exclude the specified taxa.
    if exclude_taxa:
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
    mf = mf.assign(**{'sample-id': mf.index})
    df = df.drop(columns=cols)

    df, mf = _filter_samples(df, mf, exclude_samples, include_samples)

    # If provided, only include the specified samples.
    if sample_names:
        df = df.loc[sample_names]
        mf = mf.loc[sample_names]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Add a pseudocount.
    if pseudocount:
        df = df + 1

    # Convert counts to proportions.
    df = df.div(df.sum(axis=1), axis=0)

    # Sort the columns (i.e. species) by their mean abundance.
    df = _sort_by_mean(df)

    # If provided, collapse extra species to the Others column.
    if count is not 0 and taxa_names is not None:
        m = "Cannot use 'count' and 'taxa_names' arguments together"
        raise ValueError(m)
    elif count is not 0:
        others = df.iloc[:, count-1:].sum(axis=1)
        df = df.iloc[:, :count-1]
        df['Others'] = others
    elif taxa_names is not None:
        df = df[taxa_names]
    else:
        pass

    if sort_by_names:
        df = df.reindex(sorted(df.columns), axis=1)

    _taxa_names = df.columns

    df2 = df * 100

    if hue is not None:
        df2 = pd.concat([df2, mf[hue]], axis=1, join='inner')
        df2 = pd.melt(df2, id_vars=[hue])
    else:
        df2 = pd.melt(df2)



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
        sns.swarmplot(x='variable',
                      y='value',
                      hue=hue,
                      hue_order=hue_order,
                      data=df2,
                      ax=ax,
                      color='black',
                      size=size,
                      dodge=True)
    else:
        remove_duplicates = False

    ax.set_xlabel('')
    ax.set_ylabel('Relative abundance (%)')



    # If provided, output the dataframe as a .csv file.
    if csv_file is not None:
        df3 = pd.concat([df, mf], axis=1, join='inner')
        df3.to_csv(csv_file)


    if xlabels is not None:
        xtexts = [x.get_text() for x in ax.get_xticklabels()]
        if len(xtexts) != len(xlabels):
            m = f"Expected {len(xtexts)} items, but found {len(xlabels)}"
            raise ValueError(m)
        ax.set_xticklabels(xlabels)

    a = ax.get_xticklabels()

    if brief_xlabels:
        a = [_pretty_taxa(x) for x in a]

    ax.set_xticklabels(a, rotation=45, ha='right')

    ax = _legend_handler(ax,
                         show_legend,
                         legend_loc,
                         remove_duplicates=remove_duplicates)

    ax = _artist(ax, **kwargs)

    return ax










def ancom_volcano_plot(ancom,
                       ax=None,
                       figsize=None,
                       **kwargs):
    """
    This method creates an ANCOM volcano plot.

    Parameters
    ----------
    ancom : str
        Path to the visualization file from the 'qiime composition ancom' 
        command.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

    Returns
    -------
    matplotlib.axes.Axes
        Returns the Axes object with the plot drawn onto it.
    """
    t = TemporaryDirectory()
    Visualization.load(ancom).export_data(t.name)
    df = pd.read_table(f'{t.name}/data.tsv')
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(df.clr, df.W, s=80, c='black', alpha=0.5)
    ax.set_xlabel('clr')
    ax.set_ylabel('W')

    ax = _artist(ax, **kwargs)

    return ax










# -- Other plotting methods --------------------------------------------------

def beta_2d_plot_gallery(ordination,
                         metadata,
                         prefix,
                         targets=None,
                         nrows=3,
                         ncols=4,
                         figsize=None,
                         **kwargs):

    """
    This method extends the 'beta_2d_plot' method and allows the user to 
    automatically create multiple figures of 2D beta diversity plot. This 
    method is useful when there are multiple variables to be tested.

    Parameters
    ----------
    ordination : str or qiime2.sdk.result.Artifact
        Artifact file or object from ordination.
    metadata : str or qiime2.metadata.metadata.Metadata
        Metadata file or object.
    prefix : str
        File prefix.
    targets: list, optional
        List of targeted columns. Otherwise, use all columns in the metadata.
    nrows : int, default: 3
        Number of rows of the subplot grid.
    ncols : int, default: 4
        Number of rows of the subplot grid.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    """
    if targets is None:
        mf = get_mf(metadata)
        _targets = mf.columns
    else:
        _targets = targets

    n_total = len(_targets)
    n_panels = nrows * ncols
    n_figures = math.ceil(len(_targets) / n_panels)

    i = 0 # Total number of panels.

    for j in range(n_figures):
        k = 0 # Number of panels within figure.
        filename = f'{prefix}-{j}.png'

        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize)
        
        for row in axes:
            for col in row:

                kwargs = {**kwargs,
                          'ax': col,
                          'hue': _targets[i]}

                beta_2d_plot(ordination, metadata, show_legend=True, **kwargs)

                legend_labels = [x.get_text() for x in col.legend().get_texts()] 
                if len(legend_labels) > 10:
                    col.get_legend().remove()
                
                k += 1
                i += 1

                if k == n_panels:
                    print(f"Figure saved to: {filename}")
                    plt.tight_layout()
                    plt.savefig(filename)
                    k = 0

                elif i == n_total:
                    print(f"Figure saved to: {filename}")
                    plt.tight_layout()
                    plt.savefig(filename)
                    return










def addsig(x1,
           x2,
           y,
           t='',
           h=1.0,
           lw=1.0,
           lc='black',
           tc='black'):
    """
    This method adds a signifiance annotation between two groups in a box plot.

    Parameters
    ----------
    x1 : float
        Position of the first box.
    x2 : float
        Position of the second box.
    y : float
        Bottom position of the drawing.
    t : str, default: ''
        Text.
    h : float, default: 1.0
        Height of the drawing.
    lw : float, default: 1.0
        Line width.
    lc : str, default: 'black'
        Line color.
    tc : str, default: 'black'
        Text color.
    """
    plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=lw, c=lc)
    plt.text((x1+x2)*0.5, y+h, t, ha='center', va='bottom', color=tc)









