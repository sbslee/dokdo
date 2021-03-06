import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from .common import _artist

def regplot(taxon,
            csv_file,
            subject,
            category,
            group1,
            group2,
            label=None,
            ax=None,
            figsize=None,
            artist_kwargs=None):
    """Plot relative abundance data and a linear regression model fit from
    paired samples for the given taxon.

    Parameters
    ----------
    taxon : str
        Target taxon name.
    csv_file : str
        Path to the .csv file from the `taxa_abundance_box_plot` method.
    subject : str
        Column name to indicate pair information.
    category : str
        Column name to be studied.
    group1 : str
        First group in the category column.
    group2 : str
        Second group in the category column.
    label : str
        Label to use in a legend.
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
    taxa_abundance_box_plot
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    df = pd.read_csv(csv_file)
    df = df.sort_values([subject, category])
    g1 = df[df[category] == group1]
    g2 = df[df[category] == group2]
    df = pd.DataFrame({group1: g1[taxon].to_list(),
                       group2: g2[taxon].to_list()})

    if label is None:
        _label = taxon
    else:
        _label = label

    sns.regplot(data=df, x=group1, y=group2, ax=ax, label=_label)

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {"xlabel": group1,
                     "ylabel": group2,
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
