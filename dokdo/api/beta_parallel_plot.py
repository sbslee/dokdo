from . import common

from qiime2 import Artifact

import matplotlib.pyplot as plt
import pandas as pd
from skbio.stats.ordination import OrdinationResults

def beta_parallel_plot(
    artifact, hue=None, hue_order=None, metadata=None, count=5,
    ax=None, figsize=None
):
    """
    Create a parallel plot from PCoA results.

    +---------------------+---------------------------------------------------+
    | q2-diversity plugin | Example                                           |
    +=====================+===================================================+
    | QIIME 2 CLI         | qiime diversity pcoa [OPTIONS]                    |
    +---------------------+---------------------------------------------------+
    | QIIME 2 API         | from qiime2.plugins.diversity.methods import pcoa |
    +---------------------+---------------------------------------------------+

    Parameters
    ----------
    artifact : str or qiime2.Artifact
        Artifact file or object from the q2-diversity plugin with the
        semantic type ``PCoAResults`` or
        ``PCoAResults % Properties('biplot')``.
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

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    dokdo.api.ordinate
    dokdo.api.beta_2d_plot
    dokdo.api.beta_3d_plot
    dokdo.api.beta_scree_plot

    Examples
    --------
    Below is a simple example:

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()
        qza_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/unweighted_unifrac_pcoa_results.qza'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'
        dokdo.beta_parallel_plot(qza_file,
                                 figsize=(8, 6))
        plt.tight_layout()

    .. image:: images/beta_parallel_plot-1.png

    We can group the lines by body-site:

    .. code:: python3

        dokdo.beta_parallel_plot(qza_file,
                                 metadata=metadata_file,
                                 hue='body-site',
                                 figsize=(8, 6))
        plt.tight_layout()

    .. image:: images/beta_parallel_plot-2.png
    """
    if isinstance(artifact, str):
        _pcoa_results = Artifact.load(artifact)
    else:
        _pcoa_results = artifact

    ordination_results = _pcoa_results.view(OrdinationResults)

    props = ordination_results.proportion_explained * 100
    props = [f'Axis {i+1} ({x:.2f}%)' for i, x in enumerate(props[:count])]

    df = ordination_results.samples.copy().iloc[:, :count]

    if hue is None:
        col = df.index
    else:
        mf = common.get_mf(metadata)
        col = mf[hue]

    df = df.assign(Target=col)

    if isinstance(hue_order, list):
        d = {x:i for i, x in enumerate(hue_order)}
        df = df.iloc[df['Target'].map(d).argsort()]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    pd.plotting.parallel_coordinates(df, 'Target',
        color=plt.cm.get_cmap('tab10').colors, ax=ax)

    if hue is None:
        ax.get_legend().remove()

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels(props)

    return ax
