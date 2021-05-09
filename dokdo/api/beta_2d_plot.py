from qiime2 import Artifact
import pandas as pd
import seaborn as sns
from .common import _artist
import matplotlib.pyplot as plt
from skbio.stats.ordination import OrdinationResults
import dokdo

def beta_2d_plot(
    pcoa_results, metadata=None, hue=None, size=None,
    style=None, s=80, alpha=None, ax=None,
    figsize=None, hue_order=None, style_order=None,
    legend_type='brief', artist_kwargs=None
):
    """Create a 2D scatter plot from PCoA results.

    +---------------------+---------------------------------------------------+
    | q2-diversity plugin | Example                                           |
    +=====================+===================================================+
    | QIIME 2 CLI         | qiime diversity pcoa [OPTIONS]                    |
    +---------------------+---------------------------------------------------+
    | QIIME 2 API         | from qiime2.plugins.diversity.methods import pcoa |
    +---------------------+---------------------------------------------------+

    Parameters
    ----------
    pcoa_results : str or qiime2.Artifact
        Artifact file or object corresponding to PCoAResults or
        PCoAResults % Properties('biplot').
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    hue : str, optional
        Grouping variable that will produce points with different colors.
    size : str, optional
        Grouping variable that will produce points with different sizes.
    style : str, optional
        Grouping variable that will produce points with different markers.
    s : float, default: 80.0
        Marker size.
    alpha : float, optional
        Proportional opacity of the points.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    style_order : list, optional
        Specify the order of categorical levels of the 'style' semantic.
    legend_type : str, default: 'brief'
        Legend type as in seaborn.scatterplot ('brief' or 'full').
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    ordinate
    beta_3d_plot
    beta_scree_plot
    beta_parallel_plot
    addbiplot

    Examples
    --------
    Below is a simple example.

    >>> qza_file = f'{data_dir}/moving-pictures-tutorial/unweighted_unifrac_pcoa_results.qza'
    >>> metadata_file = f'{data_dir}/moving-pictures-tutorial/sample-metadata.tsv'
    >>> dokdo.beta_2d_plot(qza_file)
    >>> plt.tight_layout()

    .. image:: images/beta_2d_plot-1.png
      :width: 600

    We can color the datapoints with ``hue``. We can also change the
    style of datapoints with ``style``. If the variable of interest is
    numeric, we can use ``size`` to control the size of datapoints.
    Finally, we can combine all those groupings.

    >>> fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, figsize=(8, 8))
    >>> artist_kwargs1 = dict(show_legend=True, title="hue='body-site'")
    >>> artist_kwargs2 = dict(show_legend=True, title="style='subject'")
    >>> artist_kwargs3 = dict(show_legend=True, title="size='days-since-experiment-start'")
    >>> artist_kwargs4 = dict(title="Combined groupings")
    >>> dokdo.beta_2d_plot(qza_file, metadata_file, ax=ax1, hue='body-site', artist_kwargs=artist_kwargs1)
    >>> dokdo.beta_2d_plot(qza_file, metadata_file, ax=ax2, style='subject', artist_kwargs=artist_kwargs2)
    >>> dokdo.beta_2d_plot(qza_file, metadata_file, ax=ax3, size='days-since-experiment-start', artist_kwargs=artist_kwargs3)
    >>> dokdo.beta_2d_plot(qza_file, metadata_file, ax=ax4, hue='body-site', style='subject', size='days-since-experiment-start', artist_kwargs=artist_kwargs4)
    >>> plt.tight_layout()

    .. image:: images/beta_2d_plot-2.png
      :width: 600
    """
    if isinstance(pcoa_results, str):
        _pcoa_results = Artifact.load(pcoa_results)
    else:
        _pcoa_results = pcoa_results

    ordination_results = _pcoa_results.view(OrdinationResults)

    df1 = ordination_results.samples.iloc[:, :2]
    df1.columns = ['A1', 'A2']

    if metadata is None:
        df2 = df1
    else:
        mf = dokdo.get_mf(metadata)
        df2 = pd.concat([df1, mf], axis=1, join='inner')

    props = ordination_results.proportion_explained

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

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': f'Axis 1 ({props[0]*100:.2f} %)',
                     'ylabel': f'Axis 2 ({props[1]*100:.2f} %)',
                     'hide_xticks': True,
                     'hide_yticks': True,
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
