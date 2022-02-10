from . import common

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from qiime2 import Artifact

def alpha_diversity_plot(
    artifact, metadata, where, add_swarmplot=False, order=None,
    hide_nsizes=False, ax=None, figsize=None
):
    """
    Create an alpha diversity plot.

    +-----------------------+--------------------------------------------------------------------------+
    | q2-diversity plugin   | Example                                                                  |
    +=======================+==========================================================================+
    | QIIME 2 CLI           | qiime diversity core-metrics-phylogenetic [OPTIONS]                      |
    +-----------------------+--------------------------------------------------------------------------+
    | QIIME 2 API           | from qiime2.plugins.diversity.pipelines import core_metrics_phylogenetic |
    +-----------------------+--------------------------------------------------------------------------+

    Parameters
    ----------
    artifact : str, qiime2.Artifact, or pandas.DataFrame
        Artifact file or object from the q2-diversity plugin with the
        semantic type ``SampleData[AlphaDiversity]``. If you are importing
        data from a software tool other than QIIME 2, then you can provide a
        :class:`pandas.DataFrame` object in which the row index is sample
        names and the only column is diversity values with its header being
        the name of metrics used (e.g. 'faith_pd').
    metadata : str or qiime2.Metadata
        Metadata file or object.
    where : str
        Column name to be used for the x-axis.
    add_swarmplot : bool, default: False
        Add a swarm plot on top of the box plot.
    order : list, optional
        Order to plot the categorical levels in.
    hide_nsizes : bool, default: False
        Hide sample size from x-axis labels.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    Examples
    --------
    Below is a simple example.

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()
        qza_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/faith_pd_vector.qza'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'
        dokdo.alpha_diversity_plot(qza_file, metadata_file, 'body-site')
        plt.tight_layout()

    .. image:: images/alpha_diversity_plot.png
    """
    if isinstance(artifact, str):
        _alpha_diversity = Artifact.load(artifact)
        df = _alpha_diversity.view(pd.Series).to_frame()
    elif isinstance(artifact, pd.DataFrame):
        df = artifact
    else:
        _alpha_diversity = artifact
        df = _alpha_diversity.view(pd.Series).to_frame()

    mf = common.get_mf(metadata)
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

    ax.set_xlabel(where)
    ax.set_ylabel(metric)

    return ax
