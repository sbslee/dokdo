import matplotlib.pyplot as plt

def addsig(x1,
           x2,
           y,
           t='',
           h=1.0,
           lw=1.0,
           lc='black',
           tc='black',
           ax=None,
           figsize=None,
           fontsize=None):
    """Add signifiance annotation between two groups in a box plot.

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
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=lw, c=lc)
    ax.text((x1+x2)*0.5, y+h, t, ha='center', va='bottom',
            color=tc, fontsize=fontsize)
    return ax
