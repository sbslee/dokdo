import tempfile
from .common import _parse_input, _artist
import pandas as pd
import skbio as sb
import matplotlib.pyplot as plt
from qiime2 import Artifact

def distance_matrix_plot(
    distance_matrix, bins=100, pairs=None, ax=None,
    figsize=None, density=False, artist_kwargs=None
):
    """Create a histogram from a distance matrix.

    Parameters
    ----------
    distance_matrix : str or qiime2.Artifact
        Artifact file or object with the semantic type
        `DistanceMatrix`.
    bins : int, optional
        Number of bins to be displayed.
    pairs : list, optional
        List of sample pairs to be shown in red vertical lines.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    density : bool, default: False
        If True, draw and return a probability density.
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    Examples
    --------
    Below is a simple example.

    >>> qza_file = f'{data_dir}/moving-pictures-tutorial/unweighted_unifrac_distance_matrix.qza'
    >>> dokdo.distance_matrix_plot(qza_file)
    >>> plt.tight_layout()

    .. image:: images/distance_matrix_plot-1.png

    We can indicate the distance between any two samples on top of the
    histogram using ``pairs``.

    >>> dokdo.distance_matrix_plot(qza_file,
    ...                            pairs=[['L1S8', 'L1S57'],
    ...                                   ['L2S175', 'L2S204']])
    >>> plt.tight_layout()

    .. image:: images/distance_matrix_plot-2.png

    Finally, we can show a density histogram.

    >>> dokdo.distance_matrix_plot(qza_file, density=True)
    >>> plt.tight_layout()

    .. image:: images/distance_matrix_plot-3.png
    """
    if isinstance(distance_matrix, str):
        _distance_matrix = Artifact.load(distance_matrix)
    else:
        _distance_matrix = distance_matrix

    dist = _distance_matrix.view(sb.DistanceMatrix)
    cdist = dist.condensed_form()

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    ax.hist(cdist, bins=bins, density=density)

    # https://stackoverflow.com/a/36867493/7481899
    def square_to_condensed(i, j, n):
        assert i != j, "no diagonal elements in condensed matrix"
        if i < j:
            i, j = j, i
        return n*j - j*(j+1)//2 + i - 1 - j

    if pairs:
        idx = []

        for pair in pairs:
            i = square_to_condensed(dist.index(pair[0]), dist.index(pair[1]), len(dist.ids))
            idx.append(cdist[i])

        for i in idx:
            ax.axvline(x=i, c='red')

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': 'Distance',
                     'ylabel': 'Frequency',
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
