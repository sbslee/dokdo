import pandas as pd
import seaborn as sns
import numpy as np

def addpairs(taxon, csv_file, subject, category, groups,
             ax=None, figsize=None, width=0.8, **kwargs):
    """Add paired lines in a plot created by `dokdo.taxa_abundance_box_plot`.

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
    groups : list
        Groups in the category column.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    width : float, default: 0.8
        Width of all the elements for one level of the grouping variable.
    kwargs : other keyword arguments
        All other keyword arguments are passed to `seaborn.lineplot`.

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
    df = df[[subject, category, taxon]]
    df = df.rename(columns={taxon: 'perc'})

    n = len(groups)
    half = width / n / 2
    positions = np.linspace(0, width - half * 2, num=n) - half * (n-1)

    d = dict(zip(groups, positions))
    df[category] = df[category].map(d)
    df = df.pivot(index=category, columns=subject, values='perc')
    ax = sns.lineplot(data=df, dashes=False, legend=False, ax=ax, **kwargs)
    ax.set_xlabel('')

    return ax
