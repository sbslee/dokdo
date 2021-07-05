import tempfile

from . import common

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def read_quality_plot(visualization, strand='forward', ax=None, figsize=None):
    """
    Create a read quality plot.

    +-----------------+--------------------------------------------------------+
    | q2-demux plugin | Example                                                |
    +=================+========================================================+
    | QIIME 2 CLI     | qiime demux summarize [OPTIONS]                        |
    +-----------------+--------------------------------------------------------+
    | QIIME 2 API     | from qiime2.plugins.demux.visualizers import summarize |
    +-----------------+--------------------------------------------------------+

    Parameters
    ----------
    visualization : str or qiime2.Visualization
        Visualization file or object from the q2-demux plugin.
    strand : {'forward', 'reverse'}, default: 'forward'
        Read strand to be displayed.
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
    Below is a simple example:

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()
        fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(10, 5))
        qzv_file = '/Users/sbslee/Desktop/dokdo/data/atacama-soil-microbiome-tutorial/demux-subsample.qzv'
        dokdo.read_quality_plot(qzv_file, strand='forward', ax=ax1)
        dokdo.read_quality_plot(qzv_file, strand='reverse', ax=ax2)
        ax1.set_title('Forward read')
        ax2.set_title('Reverse read')
        ax2.set_ylabel('')
        ax2.set_yticklabels([])
        ax2.autoscale(enable=True, axis='x', tight=False)
        plt.tight_layout()

    .. image:: images/read_quality_plot.png
    """

    with tempfile.TemporaryDirectory() as t:
        common.export(visualization, t)
        df = pd.read_table(f'{t}/{strand}-seven-number-summaries.tsv',
                           index_col=0,
                           skiprows=[1])

    df = pd.melt(df.reset_index(),
                 id_vars=['index'],
                 var_name='Base',
                 value_name='Score')
    df['Base'] = df['Base'].astype('int')

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.boxplot(x='Base',
                y='Score',
                data=df,
                ax=ax,
                fliersize=0,
                boxprops=dict(color='white', edgecolor='black'),
                medianprops=dict(color='red'),
                whiskerprops=dict(linestyle=':'))

    xticks = [int(x) for x in np.linspace(0, df['Base'].max(), 11).tolist()]

    ax.set_xlabel('Sequence base')
    ax.set_ylabel('Quality score')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks)
    ax.set_ylim([0, 45])

    return ax
