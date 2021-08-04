import matplotlib.pyplot as plt

def addsig(
    x1, x2, y, t='', h=1.0, lw=1.0, lc='black',
    tc='black', ax=None, figsize=None, fontsize=None
):
    """
    Add signifiance annotation between two groups in a box plot.

    Parameters
    ----------
    x1 : float
        Position of the first box.
    x2 : float
        Position of the second box.
    y : float
        Bottom position of the drawing.
    t : str, default: ''
        Text.
    h : float, default: 1.0
        Height of the drawing.
    lw : float, default: 1.0
        Line width.
    lc : str, default: 'black'
        Line color.
    tc : str, default: 'black'
        Text color.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    fontsize : float, optional
        Sets the fontsize.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    Examples
    --------

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()
        vector_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/faith_pd_vector.qza'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'
        ax = dokdo.alpha_diversity_plot(vector_file, metadata_file, 'body-site', figsize=(8, 5))
        dokdo.addsig(0, 1, 20, t='***', ax=ax)
        dokdo.addsig(1, 2, 26, t='ns', ax=ax)
        plt.tight_layout()

    .. image:: images/addsig.png
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=lw, c=lc)
    ax.text((x1+x2)*0.5, y+h, t, ha='center', va='bottom',
            color=tc, fontsize=fontsize)
    return ax
