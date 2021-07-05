import tempfile

from . import common

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def ancom_volcano_plot(
    visualization, ax=None, figsize=None, **kwargs
):
    """
    Create an ANCOM volcano plot.

    +-----------------------+----------------------------------------------------------+
    | q2-composition plugin | Example                                                  |
    +=======================+==========================================================+
    | QIIME 2 CLI           | qiime composition ancom [OPTIONS]                        |
    +-----------------------+----------------------------------------------------------+
    | QIIME 2 API           | from qiime2.plugins.composition.visualizers import ancom |
    +-----------------------+----------------------------------------------------------+

    Parameters
    ----------
    visualization : str or qiime2.Visualization
        Visualization file or object from the q2-composition plugin.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    kwargs
        Other keyword arguments will be passed down to
        :meth:`seaborn.scatterplot`.

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
        qzv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/ancom-subject.qzv'
        dokdo.ancom_volcano_plot(qzv_file)
        plt.tight_layout()

    .. image:: images/ancom_volcano_plot_1.png

    We can control the size, color, and transparency of data points with the
    ``s``, ``color``, and ``alpha`` options, respectively:

    .. code:: python3

        dokdo.ancom_volcano_plot(qzv_file, s=80, color='black', alpha=0.5)
        plt.tight_layout()

    .. image:: images/ancom_volcano_plot_2.png
    """
    with tempfile.TemporaryDirectory() as t:
        common.export(visualization, t)
        df = pd.read_table(f'{t}/data.tsv')

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.scatterplot(x='clr', y='W', data=df, ax=ax, **kwargs)

    ax.set_xlabel('clr')
    ax.set_ylabel('W')

    return ax
