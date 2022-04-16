from . import common

import pandas as pd
import matplotlib.pyplot as plt
from skbio.stats.ordination import OrdinationResults
from qiime2 import Artifact

def beta_3d_plot(
    artifact, metadata=None, hue=None, azim=-60, elev=30, s=80, ax=None,
    figsize=None, hue_order=None, palette=None
):
    """
    Create a 3D scatter plot from PCoA results.

    In addition to creating a PCoA plot, this method prints out the
    proportions explained by each axis.

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
        ``PCoAResults % Properties('biplot')``. If you are importing
        data from a software tool other than QIIME 2, then you can provide a
        :class:`pandas.DataFrame` object in which the row index is sample
        names and the first, second, and third columns indicate the first,
        second, and third PCoA axes, respectively.
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
    palette : dict
        Dictionary for choosing the colors to use when mapping the ``hue``
        semantic.

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
                           metadata=metadata_file,
                           hue='body-site',
                           figsize=(8, 8))
        plt.tight_layout()

    .. code-block:: text

        # Explained proportions computed by QIIME 2:
        # 33.94% for Axis 1
        # 25.90% for Axis 2
        # 6.63% for Axis 3

    .. image:: images/beta_3d_plot-1.png

    We can control the camera angle with ``elev`` and ``azim``:

    .. code:: python3

        fig = plt.figure(figsize=(14, 7))
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')
        dokdo.beta_3d_plot(qza_file,
                           metadata=metadata_file,
                           ax=ax1,
                           hue='body-site',
                           elev=15)
        dokdo.beta_3d_plot(qza_file,
                           metadata=metadata_file,
                           ax=ax2,
                           hue='body-site',
                           azim=70)
        plt.tight_layout()

    .. image:: images/beta_3d_plot-2.png

    We can control categorical mapping of the ``hue`` variable with
    ``palette``:

    .. code:: python3

        palette = {'gut': 'yellow', 'left palm': 'green', 'right palm': 'blue', 'tongue': 'red'}
        dokdo.beta_3d_plot(qza_file, metadata=metadata_file, hue='body-site', palette=palette, figsize=(8, 8))
        plt.tight_layout()

    .. image:: images/beta_3d_plot-3.png
    """
    if isinstance(artifact, pd.DataFrame):
        df = artifact
    else:
        if isinstance(artifact, str):
            _pcoa_results = Artifact.load(artifact)
        else:
            _pcoa_results = artifact
        ordination_results = _pcoa_results.view(OrdinationResults)
        df = ordination_results.samples.iloc[:, :3]
        props = ordination_results.proportion_explained
        print('# Explained proportions computed by QIIME 2:')
        print(f'# {props[0]*100:.2f}% for Axis 1')
        print(f'# {props[1]*100:.2f}% for Axis 2')
        print(f'# {props[2]*100:.2f}% for Axis 3')

    df.columns = ['Axis 1', 'Axis 2', 'Axis 3']

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
        ax.scatter(df['Axis 1'], df['Axis 2'], df['Axis 3'], s=s)
    else:
        if hue_order is None:
            _hue_order = df[hue].unique()
        else:
            _hue_order = hue_order
        for label in _hue_order:
            a = df[df[hue] == label]
            if palette is None:
                palette = {x: None for x in _hue_order}
            c = palette[label]
            ax.scatter(a['Axis 1'], a['Axis 2'], a['Axis 3'], label=label, s=s, c=c)
            ax.legend()

    ax.set_xlabel('Axis 1')
    ax.set_ylabel('Axis 2')
    ax.set_zlabel('Axis 3')

    return ax
