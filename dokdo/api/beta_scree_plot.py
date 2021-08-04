import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from skbio.stats.ordination import OrdinationResults
from qiime2 import Artifact

def beta_scree_plot(artifact, count=5, color='blue', ax=None, figsize=None):
    """
    Create a scree plot from PCoA results.

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
    count : int, default: 5
        Number of principal components to be displayed.
    color : str, default: 'blue'
        Bar color.
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
    dokdo.api.beta_parallel_plot

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
        dokdo.beta_scree_plot(qza_file)
        plt.tight_layout()

    .. image:: images/beta_scree_plot.png
    """
    if isinstance(artifact, str):
        _pcoa_results = Artifact.load(artifact)
    else:
        _pcoa_results = artifact

    ordination_results = _pcoa_results.view(OrdinationResults)
    props = ordination_results.proportion_explained
    df = pd.DataFrame({'PC': [f'Axis {x+1}' for x in range(len(props))],
                       'Proportion': props * 100})
    sliced_df = df.head(count)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.barplot(x='PC', y='Proportion', data=sliced_df, color=color, ax=ax)

    ax.set_ylabel('Variation explained (%)')

    return ax
