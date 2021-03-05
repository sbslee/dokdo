import tempfile
from ..api import _parse_input, _artist
import pandas as pd
import seaborn as sns
import numpy as np

def read_quality_plot(demux, strand='forward', ax=None,
                      figsize=None, artist_kwargs=None):
    """
    This method creates a read quality plot.

    Parameters
    ----------
    demux : str or qiime2.Visualization
        Visualization file or object from the q2-demux plugin.
    strand : str, default: 'forward'
        Read strand to be displayed (either 'forward' or 'reverse').
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    Notes
    -----
    Example usage of the q2-demux plugin:
        CLI -> $ qiime demux summarize [OPTIONS]
        API -> from qiime2.plugins.demux.visualizers import summarize
    """
    l = ['forward', 'reverse']

    if strand not in l:
        raise ValueError(f"Strand should be one of the following: {l}")

    with tempfile.TemporaryDirectory() as t:
        _parse_input(demux, t)

        df = pd.read_table(f'{t}/{strand}-seven-number-summaries.tsv',
                           index_col=0,
                           skiprows=[1])

    df = pd.melt(df.reset_index(), id_vars=['index'])
    df['variable'] = df['variable'].astype('int64')
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.boxplot(x='variable',
                y='value',
                data=df,
                ax=ax,
                fliersize=0,
                boxprops=dict(color='white', edgecolor='black'),
                medianprops=dict(color='red'),
                whiskerprops=dict(linestyle=':'))

    xticks = np.arange(df['variable'].min(), df['variable'].max(), 20).tolist()

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': 'Sequence base',
                     'ylabel': 'Quality score',
                     'xticks': xticks,
                     'xticklabels': xticks,
                     'ymin': 0,
                     'ymax': 45,
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
