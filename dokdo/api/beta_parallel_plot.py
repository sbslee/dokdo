from qiime2 import Artifact
import matplotlib.pyplot as plt
from .common import _artist
from skbio.stats.ordination import OrdinationResults
import pandas as pd
import dokdo

def beta_parallel_plot(
    pcoa_results, hue=None, hue_order=None,
    metadata=None, count=5, ax=None,
    figsize=None, artist_kwargs=None
):
    """Create a parallel plot from PCoA results.

    Parameters
    ----------
    pcoa_results : str or qiime2.Artifact
        Artifact file or object corresponding to PCoAResults.
    hue : str, optional
        Grouping variable that will produce lines with different colors.
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object. Required if 'hue' is used.
    count : int, default: 5
        Number of principal components to be displayed.
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

    See Also
    --------
    ordinate
    beta_2d_plot
    beta_3d_plot
    beta_scree_plot

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
        >>> metadata_file = f'{data_dir}/moving-pictures-tutorial/sample-metadata.tsv'
        >>> dokdo.beta_parallel_plot(qza_file)
        >>> plt.tight_layout()

    We can group the lines by body-site.

    .. plot::
        :context: close-figs

        >>> dokdo.beta_parallel_plot(qza_file,
        ...                          metadata=metadata_file,
        ...                          hue='body-site',
        ...                          artist_kwargs=dict(show_legend=True))
        >>> plt.tight_layout()
    """
    if isinstance(pcoa_results, str):
        _pcoa_results = Artifact.load(pcoa_results)
    else:
        _pcoa_results = pcoa_results

    ordination_results = _pcoa_results.view(OrdinationResults)

    props = ordination_results.proportion_explained * 100
    props = [f'Axis {i+1} ({x:.2f}%)' for i, x in enumerate(props[:count])]

    df = ordination_results.samples.copy().iloc[:, :count]

    if hue is None:
        col = df.index
    else:
        mf = dokdo.get_mf(metadata)
        col = mf[hue]

    df = df.assign(Target=col)

    if isinstance(hue_order, list):
        d = {x:i for i, x in enumerate(hue_order)}
        df = df.iloc[df['Target'].map(d).argsort()]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    pd.plotting.parallel_coordinates(df, 'Target',
        color=plt.cm.get_cmap('tab10').colors)

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': '',
                     'ylabel': '',
                     'xticklabels': props,
                     'legend_title': hue,
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
