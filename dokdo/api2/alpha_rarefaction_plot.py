import tempfile
import seaborn as sns
from ..api import _parse_input, _artist
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def alpha_rarefaction_plot(rarefaction, hue='sample-id', metric='shannon',
                           ax=None, figsize=None, hue_order=None,
                           units=None, estimator='mean', seed=1,
                           artist_kwargs=None):
    """
    This method creates an alpha rarefaction plot.

    Parameters
    ----------
    rarefaction : str or qiime2.Visualization
        Visualization file or object from the q2-diversity plugin.
    hue : str, default: 'sample-id'
        Grouping variable that will produce lines with different colors. If not
        provided, sample IDs will be used.
    metric : str, default: 'shannon'
        Diversity metric ('shannon', 'observed_features', or 'faith_pd').
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    units : str, optional
        Grouping variable identifying sampling units. When used, a separate
        line will be drawn for each unit with appropriate semantics, but no
        legend entry will be added.
    estimator : str, default: 'mean', optional
        Method for aggregating across multiple observations of the y variable
        at the same x level. If None, all observations will be drawn.
    seed : int, default: 1
        Seed for reproducible bootstrapping.
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    Notes
    -----
    Example usage of the q2-diversity plugin:
        CLI -> qiime diversity alpha-rarefaction [OPTIONS]
        API -> from qiime2.plugins.diversity.visualizers import alpha_rarefaction
    """
    l = ['observed_features', 'faith_pd', 'shannon']

    if metric not in l:
        raise ValueError(f"Metric should be one of the following: {l}")

    with tempfile.TemporaryDirectory() as t:
        _parse_input(rarefaction, t)

        df = pd.read_csv(f'{t}/{metric}.csv', index_col=0)

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
                 seed=seed)

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': 'Sequencing depth',
                     'ylabel': metric,
                     'plot_method': 'alpha_rarefaction_plot',
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
