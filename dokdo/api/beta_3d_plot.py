from . import common

import pandas as pd
import matplotlib.pyplot as plt
from skbio.stats.ordination import OrdinationResults
from qiime2 import Artifact

def beta_3d_plot(
    artifact, metadata=None, hue=None, azim=-60,
    elev=30, s=80, ax=None, figsize=None,
    hue_order=None, artist_kwargs=None
):
    """
    Create a 3D scatter plot from PCoA results.

    +---------------------+---------------------------------------------------+
    | q2-diversity plugin | Example                                           |
    +=====================+===================================================+
    | QIIME 2 CLI         | qiime diversity pcoa [OPTIONS]                    |
    +---------------------+---------------------------------------------------+
    | QIIME 2 API         | from qiime2.plugins.diversity.methods import pcoa |
    +---------------------+---------------------------------------------------+

    Parameters
    ----------
    artifact : str or qiime2.Artifact
        Artifact file or object from the q2-diversity plugin with the
        semantic type ``PCoAResults`` or
        ``PCoAResults % Properties('biplot')``.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    hue : str, optional
        Grouping variable that will produce points with different colors.
    azim : int, default: -60
        Azimuthal viewing angle.
    elev : int, default: 30
        Elevation viewing angle.
    s : float, default: 80.0
        Marker size.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    dokdo.api.ordinate
    dokdo.api.beta_2d_plot
    dokdo.api.beta_scree_plot
    dokdo.api.beta_parallel_plot
    dokdo.api.addbiplot

    Examples
    --------
    Below is a simple example:

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()
        qza_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/unweighted_unifrac_pcoa_results.qza'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'
        dokdo.beta_3d_plot(qza_file,
                           metadata_file,
                           'body-site',
                           figsize=(8, 8))
        plt.tight_layout()

    .. image:: images/beta_3d_plot-1.png

    We can control the camera angle with ``elev`` and ``azim``:

    .. code:: python3

        fig = plt.figure(figsize=(14, 7))
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        dokdo.beta_3d_plot(qza_file,
                           metadata_file,
                           ax=ax1,
                           hue='body-site',
                           elev=15)
        dokdo.beta_3d_plot(qza_file,
                           metadata_file,
                           ax=ax2,
                           hue='body-site',
                           azim=70)
        plt.tight_layout()

    .. image:: images/beta_3d_plot-2.png
    """
    if isinstance(artifact, str):
        _pcoa_results = Artifact.load(artifact)
    else:
        _pcoa_results = artifact

    ordination_results = _pcoa_results.view(OrdinationResults)

    df = ordination_results.samples.iloc[:, :3]
    df.columns = ['A1', 'A2', 'A3']

    props = ordination_results.proportion_explained

    if metadata is None:
        df = df
    else:
        mf = common.get_mf(metadata)
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

    ax.set_xlabel(f'Axis 1 ({props[0]*100:.2f} %)')
    ax.set_ylabel(f'Axis 2 ({props[1]*100:.2f} %)')
    ax.set_zlabel(f'Axis 3 ({props[2]*100:.2f} %)')

    ax.legend()

    return ax
