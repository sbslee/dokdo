from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from qiime2 import Artifact
from qiime2 import Metadata
from qiime2 import Visualization



def ancom_volcano_plot(ancom, ax=None):
    """
    This method creates an ANCOM volcano plot.

    Parameters
    ----------
    ancom : str
        Path to the visualization file from the 'qiime composition ancom' 
        command.
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
        fig, ax = plt.subplots(figsize=(15, 10))
    ax.scatter(df.clr, df.W, s=80, c='black', alpha=0.5)
    ax.set_xlabel('clr')
    ax.set_ylabel('W')



def alpha_diversity_plot(significance, where, ax=None):
    """
    This method creates an alpha diversity plot.

    Parameters
    ----------
    significance : str
        Path to the visualization file from the 'qiime diversity 
        alpha-group-significance' command.
    where : str
        Column name to be used for the x-axis.
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
        fig, ax = plt.subplots(figsize=(15, 10))
    metric = df.columns[-1]
    boxprops = dict(color='white', edgecolor='black')
    sns.boxplot(x=where, y=metric, data=df, ax=ax, boxprops=boxprops)



def read_quality_plot(demux, strand='forward', ax=None):
    """
    This method creates a read quality plot.

    Parameters
    ----------
    demux : str
        Path to the visualization file from the 'qiime demux summarize' 
        command.
    strand : str, default: 'forward'
        Read strand to be displayed (either 'forward' or 'reverse').
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
        fig, ax = plt.subplots(figsize=(15, 10))

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



def alpha_rarefaction_plot(rarefaction, where, metric='shannon', ax=None):
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
        fig, ax = plt.subplots(figsize=(15, 10))

    sns.lineplot(x='depth', y='ASV', data=data, hue=where, ax=ax,
                 err_style='bars', sort=False)

    ax.set_xlabel('Sequencing depth')
    ax.set_ylabel(metric)



def taxa_abundance_plot(taxa, level=1, by=[], ax=None):
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
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.

    Example
    -------
    api.taxa_abundance_plot('taxa-bar-plots.qzv', level=3,
                            by=['Site', 'index'])
    """
    t = TemporaryDirectory()
    Visualization.load(taxa).export_data(t.name)
    df = pd.read_csv(f'{t.name}/level-{level}.csv', index_col=0)

    if by:
        df = df.sort_values(by=by)

    dropped = []
    for column in df.columns:
        if 'Unassigned' in column:
            continue
        elif '__' in column:
            continue
        else:
            dropped.append(column)

    df = df.drop(columns=dropped)
    df = df.T
    df = df / df.sum()
    df = df.T
    df = df.loc[:, df.mean().sort_values(ascending=False).index]

    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 10))

    df.plot.bar(stacked=True, legend=False, ax=ax,
                color=plt.cm.get_cmap('Accent').colors)

    ax.set_xlabel('Samples')
    ax.set_ylabel('Relative frequency')



def beta_2d_plot(ordination, metadata, where, ax=None):
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
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
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
        fig, ax = plt.subplots(figsize=(15, 15))

    ax.set_xlabel(f'Axis 1 ({v[0]} %)')
    ax.set_ylabel(f'Axis 2 ({v[1]} %)')
    ax.set_xticks([])
    ax.set_yticks([])

    for c in sorted(df2[where].unique()):
        i = df2[where] == c
        ax.scatter(df1[i].iloc[:, 0], df1[i].iloc[:, 1], label=c, s=80)



def beta_3d_plot(ordination, metadata, where, ax=None, azim=-60, elev=30):
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
    ax : matplotlib Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    azim : int, default: -60
        Elevation viewing angle.
    elev : int, default: 30
        Azimuthal viewing angle.
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
        fig = plt.figure(figsize=(15, 15))
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
                   df1[i].iloc[:, 2], label=c, s=80)