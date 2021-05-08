import tempfile
import pandas as pd
from .common import _parse_input, _artist
import seaborn as sns
import matplotlib.pyplot as plt

def ancom_volcano_plot(
    ancom, ax=None, figsize=None,
    s=80, artist_kwargs=None
):
    """Create an ANCOM volcano plot.

    +-----------------------+----------------------------------------------------------+
    | q2-composition plugin | Example                                                  |
    +=======================+==========================================================+
    | QIIME 2 CLI           | qiime composition ancom [OPTIONS]                        |
    +-----------------------+----------------------------------------------------------+
    | QIIME 2 API           | from qiime2.plugins.composition.visualizers import ancom |
    +-----------------------+----------------------------------------------------------+

    Parameters
    ----------
    ancom : str
        Visualization file or object from the q2-composition plugin.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    s : float, default: 80.0
        Marker size.
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    Examples
    --------
    Below is a simple example.

    .. plot::
        :context: close-figs

        >>> qzv_file = f'{data_dir}/moving-pictures-tutorial/ancom-subject.qzv'
        >>> dokdo.ancom_volcano_plot(qzv_file, figsize=(8, 5))
        >>> plt.tight_layout()
    """
    with tempfile.TemporaryDirectory() as t:
        _parse_input(ancom, t)
        df = pd.read_table(f'{t}/data.tsv')

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(data=df, x='clr', y='W', ax=ax, s=s, alpha=0.5,
                    color='black')

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': 'clr',
                     'ylabel': 'W',
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
