import tempfile

from . import common

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def denoising_stats_plot(
    artifact, metadata, where, pseudocount=False, order=None,
    hide_nsizes=False, ax=None, figsize=None,
):
    """
    Create a grouped box plot for denoising statistics from DADA2.

    +-----------------+---------------------------------------------------------+
    | q2-dada2 plugin | Example                                                 |
    +=================+=========================================================+
    | QIIME 2 CLI     | qiime dada2 denoise-paired [OPTIONS]                    |
    +-----------------+---------------------------------------------------------+
    | QIIME 2 API     | from qiime2.plugins.dada2.methods import denoise_paired |
    +-----------------+---------------------------------------------------------+

    Parameters
    ----------
    artifact : str or qiime2.Artifact
        Artifact file or object from the q2-dada2 plugin.
    metadata : str or qiime2.Metadata
        Metadata file or object.
    where : str
        Column name of the sample metadata.
    pseudocount : bool, default: False
        Add pseudocount to remove zeros.
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
        qza_file = '/Users/sbslee/Desktop/dokdo/data/atacama-soil-microbiome-tutorial/denoising-stats.qza'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/atacama-soil-microbiome-tutorial/sample-metadata.tsv'
        dokdo.denoising_stats_plot(qza_file,
                                   metadata_file,
                                   'transect-name',
                                   figsize=(8, 6))
        plt.tight_layout()

    .. image:: images/denoising_stats_plot.png
    """
    with tempfile.TemporaryDirectory() as t:
        common.export(artifact, t)
        df1 = pd.read_table(f'{t}/stats.tsv', skiprows=[1], index_col=0)

    mf = common.get_mf(metadata)

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

    ax.set_xlabel(where)
    ax.set_ylabel('Read depth')

    return ax
