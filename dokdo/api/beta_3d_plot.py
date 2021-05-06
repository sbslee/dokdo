from qiime2 import Artifact
import pandas as pd
from .common import _artist
import dokdo
import matplotlib.pyplot as plt
from skbio.stats.ordination import OrdinationResults

def beta_3d_plot(pcoa_results,
                 metadata=None,
                 hue=None,
                 azim=-60,
                 elev=30,
                 s=80,
                 ax=None,
                 figsize=None,
                 hue_order=None,
                 artist_kwargs=None):
    """Create a 3D scatter plot from PCoA results.

    Parameters
    ----------
    pcoa_results : str or qiime2.Artifact
        Artifact file or object corresponding to PCoAResults or
        PCoAResults % Properties('biplot').
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    hue : str, optional
        Grouping variable that will produce points with different colors.
    azim : int, default: -60
        Elevation viewing angle.
    elev : int, default: 30
        Azimuthal viewing angle.
    s : float, default: 80.0
        Marker size.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    ordinate
    beta_2d_plot
    beta_scree_plot
    beta_parallel_plot
    addbiplot

    Notes
    -----
    Example usage of the q2-diversity plugin:
        CLI -> qiime diversity pcoa [OPTIONS]
        API -> from qiime2.plugins.diversity.methods import pcoa

    Examples
    --------
    Below is a simple example.

    .. plot::
        :context: close-figs

        >>> import dokdo
        >>> import seaborn as sns
        >>> import matplotlib.pyplot as plt
        >>> sns.set()
        >>> data_dir = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial'
        >>> qza_file = f'{data_dir}/unweighted_unifrac_pcoa_results.qza'
        >>> metadata_file = f'{data_dir}/sample-metadata.tsv'
        >>> dokdo.beta_3d_plot(qza_file,
        ...                    metadata_file,
        ...                    'body-site',
        ...                    figsize=(6, 6),
        ...                    artist_kwargs=dict(show_legend=True))
        >>> plt.tight_layout()

    We can control the camera angle with ``elev`` and ``azim``.

    .. plot::
        :context: close-figs

        >>> fig = plt.figure(figsize=(12, 6))
        >>> ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        >>> ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        >>> dokdo.beta_3d_plot(qza_file, metadata_file, ax=ax1, hue='body-site', elev=15)
        >>> dokdo.beta_3d_plot(qza_file, metadata_file, ax=ax2, hue='body-site', azim=70)
        >>> plt.tight_layout()
    """
    if isinstance(pcoa_results, str):
        _pcoa_results = Artifact.load(pcoa_results)
    else:
        _pcoa_results = pcoa_results

    ordination_results = _pcoa_results.view(OrdinationResults)

    df = ordination_results.samples.iloc[:, :3]
    df.columns = ['A1', 'A2', 'A3']

    props = ordination_results.proportion_explained

    if metadata is None:
        df = df
    else:
        mf = dokdo.get_mf(metadata)
        df = pd.concat([df, mf], axis=1, join='inner')

    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection='3d')

    ax.view_init(azim=azim, elev=elev)

    if hue is None:
        ax.scatter(df['A1'], df['A2'], df['A3'], s=s)
    else:
        if hue_order is None:
            _hue_order = df[hue].unique()
        else:
            _hue_order = hue_order
        for label in _hue_order:
            a = df[df[hue] == label]
            ax.scatter(a['A1'], a['A2'], a['A3'], label=label, s=s)

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {'xlabel': f'Axis 1 ({props[0]*100:.2f} %)',
                     'ylabel': f'Axis 2 ({props[1]*100:.2f} %)',
                     'zlabel': f'Axis 3 ({props[2]*100:.2f} %)',
                     'hide_xticks': True,
                     'hide_yticks': True,
                     'hide_zticks': True,
                     'legend_title': hue,
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
