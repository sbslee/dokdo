from qiime2 import Artifact
from skbio.stats.ordination import OrdinationResults
import seaborn as sns
import pandas as pd
from .common import _artist
import matplotlib.pyplot as plt

def beta_scree_plot(
    pcoa_results, count=5, ax=None,
    figsize=None, color='blue', artist_kwargs=None
):
    """Create a scree plot from PCoA results.

    Parameters
    ----------
    pcoa_results : str or qiime2.Artifact
        Artifact file or object corresponding to PCoAResults.
    count : int, default: 5
        Number of principal components to be displayed.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    color : str, default: 'blue'
        Bar color.
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    ordinate
    beta_2d_plot
    beta_3d_plot
    beta_parallel_plot

    Notes
    -----
    Example usage of the q2-diversity plugin:
        CLI -> qiime diversity pcoa [OPTIONS]
        API -> from qiime2.plugins.diversity.methods import pcoa

    Examples
    --------
    Below is a simple example.

    .. plot::
        :context: close-figs

        >>> qza_file = f'{data_dir}/moving-pictures-tutorial/unweighted_unifrac_pcoa_results.qza'
        >>> dokdo.beta_scree_plot(qza_file)
        >>> plt.tight_layout()
    """
    if isinstance(pcoa_results, str):
        _pcoa_results = Artifact.load(pcoa_results)
    else:
        _pcoa_results = pcoa_results

    ordination_results = _pcoa_results.view(OrdinationResults)
    props = ordination_results.proportion_explained
    df = pd.DataFrame({'PC': [f'Axis {x+1}' for x in range(len(props))],
                       'Proportion': props * 100})
    sliced_df = df.head(count)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.barplot(x='PC', y='Proportion', data=sliced_df, color=color, ax=ax)

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': '',
                     'ylabel': 'Variation explained (%)',
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
