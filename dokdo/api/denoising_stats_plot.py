import tempfile
from .common import _parse_input, _artist
import dokdo
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def denoising_stats_plot(stats, metadata, where, ax=None, figsize=None,
                         pseudocount=False, order=None, hide_nsizes=False,
                         artist_kwargs=None):
    """Create a grouped box chart for denoising statistics from DADA 2.

    Parameters
    ----------
    stats : str or qiime2.Artifact
        Artifact file or object from the q2-dada2 plugin.
    metadata : str or qiime2.Metadata
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
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    Notes
    -----
    Example usage of the q2-dada2 plugin:
        CLI -> qiime dada2 denoise-paired [OPTIONS]
        API -> from qiime2.plugins.dada2.methods import denoise_paired
    """
    with tempfile.TemporaryDirectory() as t:
        _parse_input(stats, t)

        df1 = pd.read_table(f'{t}/stats.tsv', skiprows=[1], index_col=0)

    mf = dokdo.get_mf(metadata)

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

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': where,
                     'ylabel': 'Read depth',
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
