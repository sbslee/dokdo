from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import skbio as sb

from qiime2 import Artifact
from qiime2 import Metadata
from qiime2 import Visualization



def ancom_volcano_plot(ancom, figsize=None, ax=None):
    """
    This method creates an ANCOM volcano plot.

    Parameters
    ----------
    ancom : str
        Path to the visualization file from the 'qiime composition ancom' 
        command.
    figsize : tuple of float, optional
        Width, height in inches.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.

    Example
    -------
    api.ancom_volcano_plot('ancom-Site.qzv')
    """
    t = TemporaryDirectory()
    Visualization.load(ancom).export_data(t.name)
    df = pd.read_table(f'{t.name}/data.tsv')
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(df.clr, df.W, s=80, c='black', alpha=0.5)
    ax.set_xlabel('clr')
    ax.set_ylabel('W')



def alpha_diversity_plot(significance, where, figsize=None, ax=None):
    """
    This method creates an alpha diversity plot.

    Parameters
    ----------
    significance : str
        Path to the visualization file from the 'qiime diversity 
        alpha-group-significance' command.
    where : str
        Column name to be used for the x-axis.
    figsize : tuple of float, optional
        Width, height in inches.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.

    Example
    -------
    api.alpha_diversity_plot('shannon_group-significance.qzv', 'Site')
    """
    t = TemporaryDirectory()
    Visualization.load(significance).export_data(t.name)
    df = Metadata.load(f'{t.name}/metadata.tsv').to_dataframe()
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    metric = df.columns[-1]
    boxprops = dict(color='white', edgecolor='black')
    sns.boxplot(x=where, y=metric, data=df, ax=ax, boxprops=boxprops)



def read_quality_plot(demux, strand='forward', figsize=None, ax=None):
    """
    This method creates a read quality plot.

    Parameters
    ----------
    demux : str
        Path to the visualization file from the 'qiime demux summarize' 
        command.
    strand : str, default: 'forward'
        Read strand to be displayed (either 'forward' or 'reverse').
    figsize : tuple of float, optional
        Width, height in inches.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.

    Example
    -------
    api.read_quality_plot('demux.qzv', 'reverse')
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



def alpha_rarefaction_plot(rarefaction, where, metric='shannon',
                           figsize=None, ax=None):
    """
    This method creates an alpha rarefaction plot.

    Parameters
    ----------
    rarefaction : str
        Path to the visualization file from the 'qiime diversity 
        alpha-rarefaction' command.
    where : str
        Column name of the sample metadata.
    metric : str, default: 'shannon'
        Desired diversity metric to be displayed (either 'observed_features', 
        'faith_pd' or 'shannon').
    figsize : tuple of float, optional
        Width, height in inches.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.

    Example
    -------
    api.alpha_rarefaction_plot('alpha-rarefaction.qzv', 'Site', 
                               metric='observed_features')
    """
    t = TemporaryDirectory()
    Visualization.load(rarefaction).export_data(t.name)

    l = ['observed_features', 'faith_pd', 'shannon']
    if metric not in l:
        raise ValueError(f"Metric should be one of the following: {l}")

    df = pd.read_csv(f'{t.name}/{metric}.csv', index_col=0)
    cols = [x for x in df.columns if 'iter' not in x]
    mean = df[cols]
    data = pd.DataFrame(columns=cols)
    depths = [col.split('_')[0] for col in df.columns if 'depth' in col]
    df = df.drop(cols, axis=1)
    df.columns = depths

    for depth in depths:
        mean['ASV'] = df[depth].mean(axis=1)
        mean['depth']= depth.split('-')[-1]
        data = pd.concat([data, mean], sort=True)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.lineplot(x='depth', y='ASV', data=data, hue=where, ax=ax,
                 err_style='bars', sort=False)

    ax.set_xlabel('Sequencing depth')
    ax.set_ylabel(metric)



def taxa_abundance_plot(taxa, level=1, by=[], figsize=None, ax=None, 
                        width=0.8, count=0,
                        exclude_samples={}, exclude_taxa=[],
                        show_legend=False,
                        legend_short=False,
                        legend_loc='best', csv_file=None,
                        sort_by_names=False,
                        colors=[],
                        hide_labels=False,
                        label_columns=[]):
    """
    This method creates a taxa abundance plot.

    Parameters
    ----------
    taxa : str
        Path to the visualization file from the 'qiime taxa barplot'.
    level : int
        Taxonomic level at which the features should be collapsed.
    by : list of str
        Column name(s) to be used for sorting the samples. Using 'index' will 
        sort the samples by their name, in addition to other column name(s) 
        that may have been provided. If multiple items are provided, sorting 
        will occur by the order of the items.
    figsize : tuple of float, optional
        Width, height in inches.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    width : float
        The width of the bars.
    count : int, default: 0
        The number of taxa to display. When 0, display all.
    exclude_samples : dict
        Dictionary of column name(s) to list(s) of column value(s) to use to 
        exclude samples.
    exclude_taxa : list
        The taxa names to be excluded when matched. Case insenstivie.
    show_legend : bool, default: False
        Show the legend.
    legend_short : bool
        If true, only display the smallest taxa rank in the legend.
    legend_loc : str
        The location of the legend. Valid options include 'upper right', 
        'upper left', etc. See the 'matplotlib.pyplot.legend' method for 
        the complete list of options.
    csv_file : str, optional
        The path to the csv file.
    sort_by_names : bool
        If true, sort the columns (i.e. species) to be displayed by name.
    colors : list
        The bar colors.
    hide_labels : bool, default: False
        Hide all the x-axis labels.
    label_columns : list
        The column names to be used as the x-axis labels.
    """
    t = TemporaryDirectory()
    Visualization.load(taxa).export_data(t.name)
    df = pd.read_csv(f'{t.name}/level-{level}.csv', index_col=0)

    # If provided, sort the samples for display in the x-axis.
    if by:
        df = df.sort_values(by=by)

    # If provided, exclude the specified samples.
    if exclude_samples:
        for x in exclude_samples:
            for y in exclude_samples[x]:
                df = df[df[x] != y]

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
    dropped = []
    for column in df.columns:
        if 'Unassigned' in column:
            continue
        elif '__' in column:
            continue
        else:
            dropped.append(column)

    mf = df[dropped]
    df = df.drop(columns=dropped)




    # Convert counts to proportions.
    df = df.div(df.sum(axis=1), axis=0)

    # Sort the columns (i.e. species) by their mean abundance.
    df = df.loc[:, df.mean().sort_values(ascending=False).index]

    # If provided, collapse extra species to the Other column.
    if count is not 0:
        other = df.iloc[:, count-1:].sum(axis=1)
        df = df.iloc[:, :count-1]
        df['Other'] = other

    if sort_by_names:
        df = df.reindex(sorted(df.columns), axis=1)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if show_legend and legend_short:
        def f(s):
            ranks = s.split(';')
            for rank in reversed(ranks):
                if rank != '__':
                    x = rank
                    break
            return x
        df.columns = [f(x) for x in df.columns]

    if colors:
        c = colors
    else:
        c = plt.cm.get_cmap('Accent').colors


    df.plot.bar(stacked=True, legend=False, ax=ax, width=width, color=c)

    ax.set_xlabel('')
    ax.set_ylabel('Relative frequency')


    # Control the x-axis labels.
    if hide_labels:
        ax.set_xticks([])
    elif label_columns:
        f = lambda row: ' : '.join(row.values.astype(str))
        new_labels = mf[label_columns].apply(f, axis=1)
        ax.set_xticklabels(new_labels)
    else:
        pass


    # Control the legend.
    if show_legend:
        ax.legend(loc=legend_loc)


    if csv_file:
        df.to_csv(csv_file)





def beta_2d_plot(ordination, metadata, where, s=80, remove=[], small=[],
                 ax=None, figsize=None, show_legend=False, legend_loc='best'):
    """
    This method creates a 2D beta diversity plot.

    Parameters
    ----------
    ordination : str
        Path to the artifact file from ordination (e.g. 
        bray_curtis_pcoa_results.qza).
    metadata : str
        Path to the sample-metadata.tsv file.
    where : str
        Column name of the sample metadata.
    s : int, default: 80
        Marker size.
    remove : list of str
        Values in the column which should not be drawn when matached.
    small : dict of str
        Values in the column which should be drawn smaller when matached.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple of float, optional
        Width, height in inches.
    legend_short : bool
        If true, only display the smallest taxa rank in the legend.
    legend_loc : str
        The location of the legend. Valid options include 'upper right', 
        'upper left', etc. See the 'matplotlib.pyplot.legend' method for 
        the complete list of options.
    """
    t = TemporaryDirectory()
    Artifact.load(ordination).export_data(t.name)

    df1 = pd.read_table(f'{t.name}/ordination.txt', header=None, index_col=0,
                        skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                        skipfooter=4, engine='python', usecols=[0, 1, 2])

    df2 = Metadata.load(metadata).to_dataframe()

    df3 = pd.concat([df1, df2], axis=1, join='inner')

    f = open(f'{t.name}/ordination.txt')
    v = [round(float(x) * 100, 2) for x in f.readlines()[4].split('\t')]
    f.close()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    ax.set_xlabel(f'Axis 1 ({v[0]} %)')
    ax.set_ylabel(f'Axis 2 ({v[1]} %)')
    ax.set_xticks([])
    ax.set_yticks([])

    for c in sorted(df3[where].unique()):
        if c in remove:
            continue
        i = df3[where] == c

        _s = s / 10 if c in small else s

        ax.scatter(df3[i].iloc[:, 0], df3[i].iloc[:, 1], label=c, s=_s)

    # Control the legend.
    if show_legend:
        ax.legend(loc=legend_loc)



def beta_3d_plot(ordination, metadata, where, azim=-60, elev=30, s=80, 
                 figsize=None, ax=None, show_legend=False, legend_loc='best'):
    """
    This method creates a 3D beta diversity plot.

    Parameters
    ----------
    ordination : str
        Path to the artifact file from ordination (e.g. 
        bray_curtis_pcoa_results.qza).
    metadata : str
        Path to the sample-metadata.tsv file.
    where : str
        Column name of the sample metadata.
    azim : int, default: -60
        Elevation viewing angle.
    elev : int, default: 30
        Azimuthal viewing angle.
    s : int, default: 80
        Marker size.
    figsize : tuple of float, optional
        Width, height in inches.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    legend_short : bool
        If true, only display the smallest taxa rank in the legend.
    legend_loc : str
        The location of the legend. Valid options include 'upper right', 
        'upper left', etc. See the 'matplotlib.pyplot.legend' method for 
        the complete list of options.
    """

    t = TemporaryDirectory()
    Artifact.load(ordination).export_data(t.name)

    df1 = pd.read_table(f'{t.name}/ordination.txt', header=None, index_col=0,
                        skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                        skipfooter=4, engine='python')
    df1 = df1.sort_index()

    df2 = Metadata.load(metadata).to_dataframe()
    df2 = df2.sort_index()

    f = open(f'{t.name}/ordination.txt')
    v = [round(float(x) * 100, 2) for x in f.readlines()[4].split('\t')]
    f.close()

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

    for c in sorted(df2[where].unique()):
        i = df2[where] == c
        ax.scatter(df1[i].iloc[:, 0], df1[i].iloc[:, 1],
                   df1[i].iloc[:, 2], label=c, s=s)

    # Control the legend.
    if show_legend:
        ax.legend(loc=legend_loc)



def distance_matrix_plot(distance_matrix, bins=100, pairs={}, figsize=None, 
                         ax=None):
    """
    This method creates a histogram from a distance matrix.

    Parameters
    ----------
    distance_matrix : str
         Path to the artifact file from distance matrix computation. For 
         example, it could be an artifact from the 'qiime diversity-lib 
         jaccard' command.
    bins : int, optional
        Number of bins to be displayed.
    pairs : dict of str to list of str, optional
        Dictionary of sample pairs to be shown in red vertical lines. Keys 
        do not matter, but values have to be a list of two sample IDs.
    figsize : tuple of float, optional
        Width, height in inches.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.

    Example
    -------
    api.distance_matrix_plot('distance_matrix.qza')
    """
    t = TemporaryDirectory()
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

        for pairid, l in pairs.items():
            i = square_to_condensed(dist.index(l[0]), dist.index(l[1]), len(dist.ids))
            idx.append(cdist[i])

        for i in idx:
            ax.axvline(x=i, c='red')



def denoising_stats_plot(stats, metadata, where, figsize=None, ax=None):
    """
    This method creates a grouped box plot using denoising statistics from 
    DADA2 (i.e. the 'qiime dada2 denoise-paired' command).

    Parameters
    ----------
    stats : str
        Path to the denoising-stats.qza file.
    metadata : str
        Path to the sample-metadata.tsv file.
    where : str
        Column name of the sample metadata.
    figsize : tuple of float, optional
        Width, height in inches.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.

    Example
    -------
    api.denoising_stats_plot("denoising-stats.qza", "sample-metadata.tsv",
                             "Site")
    """
    t = TemporaryDirectory()
    Artifact.load(stats).export_data(t.name)
    df1 = pd.read_table(f'{t.name}/stats.tsv', skiprows=[1], index_col=0)
    df2 = Metadata.load(metadata).to_dataframe()
    df3 = pd.concat([df1, df2], axis=1, join='inner')
    dict = df3[where].value_counts().to_dict()
    for k, v in dict.items():
        dict[k] = f"{k} ({v})"
    df3[where].replace(dict, inplace=True)
    a = ['input', 'filtered', 'denoised', 'merged', 'non-chimeric', where]
    df4 = pd.melt(df3[a], id_vars=[where])
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(x=where, y='value', data=df4, hue='variable', ax=ax)



def paired_abundance_plot(taxa, x, y, hue, level=1, by=[],
                          exclude_samples={}, ax=None, figsize=None):
    """
    This method creates a line plot showing relative frequency of a specific 
    taxon for paired samples.

    Parameters
    ----------
    taxa : str
        Path to the visualization file from the 'qiime taxa barplot'.
    x : str
        The column to be used for the x-axis.
    y : str
        The column to be used for the y-axis (i.e. the taxon name).
    hue : str
        The column to be used for distinguishing the individual samples
    level : int
        Taxonomic level at which the features should be collapsed.
        within a pair.
    by : list of str
        Column name(s) to be used for sorting the samples. Using 'index' will 
        sort the samples by their name, in addition to other column name(s) 
        that may have been provided. If multiple items are provided, sorting 
        will occur by the order of the items.
    exclude_samples : dict of str to list of str
        Dictionary of column name(s) to list(s) of column value(s) to use to 
        exclude samples.
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple of float, optional
        Width, height in inches.
    """
    t = TemporaryDirectory()
    Visualization.load(taxa).export_data(t.name)
    df = pd.read_csv(f'{t.name}/level-{level}.csv', index_col=0)

    # If provided, sort the samples for display in the x-axis.
    if by:
        df = df.sort_values(by=[hue] + by)

    # Remove the metadata columns.
    dropped = []
    for column in df.columns:
        if 'Unassigned' in column:
            continue
        elif '__' in column:
            continue
        else:
            dropped.append(column)
    mf = df[dropped]
    df = df.drop(columns=dropped)

    # Convert counts to proportions.
    df = df.div(df.sum(axis=1), axis=0)
    mf[y] = df[y]

    # If provided, exclude the specified samples.
    if exclude_samples:
        for xx in exclude_samples:
            for yy in exclude_samples[xx]:
                mf = mf[mf[xx] != yy]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    mf.rename(columns={y: 'Y'}, inplace=True)
    mf = mf[[x, 'Y', hue] + by]

    if by:
        temp = mf[by].apply(lambda row: ', '.join(row.values.astype(str)), axis=1)

        mf[x] = mf[x] + ' (' + temp + ')'

    sns.lineplot(data=mf, x=x, y='Y', hue=hue, marker='o', sort=False, ax=ax)
    ax.tick_params(axis='x', labelrotation=90)

    ax.set_xlabel('')
    ax.set_ylabel('Relative frequency')


