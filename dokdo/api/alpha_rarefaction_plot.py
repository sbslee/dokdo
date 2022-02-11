import tempfile

from . import common

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def alpha_rarefaction_plot(
    visualization, hue='sample-id', metric='shannon', hue_order=None,
    units=None, estimator='mean', legend='brief', ax=None, figsize=None
):
    """
    Create an alpha rarefaction plot.

    +-----------------------+--------------------------------------------------------------------+
    | q2-diversity plugin   | Example                                                            |
    +=======================+====================================================================+
    | QIIME 2 CLI           | qiime diversity alpha-rarefaction [OPTIONS]                        |
    +-----------------------+--------------------------------------------------------------------+
    | QIIME 2 API           | from qiime2.plugins.diversity.visualizers import alpha_rarefaction |
    +-----------------------+--------------------------------------------------------------------+

    Parameters
    ----------
    visualization : str or qiime2.Visualization
        Visualization file or object from the q2-diversity plugin.
    hue : str, default: 'sample-id'
        Grouping variable that will produce lines with different colors. If not
        provided, sample IDs will be used.
    metric : str, default: 'shannon'
        Diversity metric ('shannon', 'observed_features', or 'faith_pd').
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    units : str, optional
        Grouping variable identifying sampling units. When used, a separate
        line will be drawn for each unit with appropriate semantics, but no
        legend entry will be added.
    estimator : str, default: 'mean', optional
        Method for aggregating across multiple observations of the y variable
        at the same x level. If None, all observations will be drawn.
    legend : str, default: 'brief'
        Legend type as in :meth:`seaborn.lineplot`.
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
        qzv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/alpha-rarefaction.qzv'
        ax = dokdo.alpha_rarefaction_plot(qzv_file,
                                          figsize=(9, 6))
        ax.legend(ncol=5)
        plt.tight_layout()

    .. image:: images/alpha_rarefaction_plot-1.png

    We can group the samples by body-site:

    .. code:: python3

        dokdo.alpha_rarefaction_plot(qzv_file,
                                     hue='body-site',
                                     metric='observed_features',
                                     figsize=(9, 6),
                                     units='sample-id',
                                     estimator=None)
        plt.tight_layout()

    .. image:: images/alpha_rarefaction_plot-2.png

    Alternatively, we can aggregate the samples by body-site:

    .. code:: python3

        dokdo.alpha_rarefaction_plot(qzv_file,
                                     hue='body-site',
                                     metric='observed_features',
                                     figsize=(9, 6))
        plt.tight_layout()

    .. image:: images/alpha_rarefaction_plot-3.png
    """
    l = ['observed_features', 'faith_pd', 'shannon']

    if metric not in l:
        raise ValueError(f"Metric should be one of the following: {l}")

    with tempfile.TemporaryDirectory() as t:
        common.export(visualization, t)
        df = pd.read_csv(f'{t}/{metric}.csv', index_col=0,
            keep_default_na=False, na_values=[''])
        df.to_csv('test.csv')

    metadata_columns = [x for x in df.columns if 'iter' not in x]

    df = pd.melt(df.reset_index(), id_vars=['sample-id'] + metadata_columns)

    df['variable'] = df['variable'].str.split('_').str[0].str.replace(
                         'depth-', '').astype(int)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.lineplot(x='variable',
                 y='value',
                 data=df,
                 hue=hue,
                 ax=ax,
                 err_style='bars',
                 sort=False,
                 units=units,
                 estimator=estimator,
                 hue_order=hue_order,
                 legend=legend)

    ax.set_xlabel('Sequencing depth')
    ax.set_ylabel(metric)

    return ax
