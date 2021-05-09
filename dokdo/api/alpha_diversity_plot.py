import pandas as pd
import seaborn as sns
from qiime2 import Artifact
import matplotlib.pyplot as plt
from .common import _artist
import dokdo

def alpha_diversity_plot(
    alpha_diversity, metadata, where,
    ax=None, figsize=None, add_swarmplot=False,
    order=None, hide_nsizes=False, artist_kwargs=None
):
    """Create an alpha diversity plot.

    Parameters
    ----------
    alpha_diversity : str or qiime2.Artifact
        Artifact file or object with the semantic type
        `SampleData[AlphaDiversity]`.
    metadata : str or qiime2.Metadata
        Metadata file or object.
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
    hide_nsizes : bool, default: False
        Hide sample size from x-axis labels.
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    Examples
    --------
    Below is a simple example.

    >>> qzv_file = f'{data_dir}/moving-pictures-tutorial/faith_pd_vector.qza'
    >>> metadata_file = f'{data_dir}/moving-pictures-tutorial/sample-metadata.tsv'
    >>> dokdo.alpha_diversity_plot(qzv_file, metadata_file, 'body-site')
    >>> plt.tight_layout()

    .. image:: images/alpha_diversity_plot.png
    """
    if isinstance(alpha_diversity, str):
        _alpha_diversity = Artifact.load(alpha_diversity)
    else:
        _alpha_diversity = alpha_diversity

    df = _alpha_diversity.view(pd.Series).to_frame()

    mf = dokdo.get_mf(metadata)
    df = pd.concat([df, mf], axis=1, join='inner')

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    metric = df.columns[0]

    boxprops = dict(color='white', edgecolor='black')

    d = {'x': where, 'y': metric, 'ax': ax, 'order': order, 'data': df}

    sns.boxplot(boxprops=boxprops, **d)

    if add_swarmplot:
        sns.swarmplot(**d)

    if hide_nsizes is False:
        nsizes = df[where].value_counts().to_dict()
        xtexts = [x.get_text() for x in ax.get_xticklabels()]
        xtexts = [f'{x} ({nsizes[x]})' for x in xtexts]
        ax.set_xticklabels(xtexts)

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': where,
                     'ylabel': metric,
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
