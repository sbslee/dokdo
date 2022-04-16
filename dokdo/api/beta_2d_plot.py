from . import common

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from skbio.stats.ordination import OrdinationResults
from qiime2 import Artifact

def beta_2d_plot(
    artifact, metadata=None, hue=None, size=None,
    style=None, s=80, alpha=None, hue_order=None, style_order=None,
    legend='brief', ax=None, figsize=None, palette=None, **kwargs
):
    """
    Create a 2D scatter plot from PCoA results.

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
    artifact : str, qiime2.Artifact, or pandas.DataFrame
        Artifact file or object from the q2-diversity plugin with the
        semantic type ``PCoAResults`` or
        ``PCoAResults % Properties('biplot')``. If you are importing
        data from a software tool other than QIIME 2, then you can provide a
        :class:`pandas.DataFrame` object in which the row index is sample
        names and the first and second columns indicate the first and second
        PCoA axes, respectively.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    hue : str, optional
        Grouping variable that will produce points with different colors.
    size : str, optional
        Grouping variable that will produce points with different sizes.
    style : str, optional
        Grouping variable that will produce points with different markers.
    s : float, default: 80.0
        Marker size.
    alpha : float, optional
        Proportional opacity of the points.
    hue_order : list, optional
        Specify the order of categorical levels of the 'hue' semantic.
    style_order : list, optional
        Specify the order of categorical levels of the 'style' semantic.
    legend : str, default: 'brief'
        Legend type as in :meth:`seaborn.scatterplot`.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    palette : string, list, dict, or matplotlib.colors.Colormap
        Method for choosing the colors to use when mapping the ``hue``
        semantic. List or dict values imply categorical mapping, while a
        colormap object implies numeric mapping.
    kwargs : other keyword arguments
        Other keyword arguments will be passed down to
        :meth:`seaborn.scatterplot()`.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    dokdo.api.ordinate
    dokdo.api.beta_3d_plot
    dokdo.api.beta_scree_plot
    dokdo.api.beta_parallel_plot
    dokdo.api.addbiplot

    Examples
    --------
    Below is a simple example.

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()
        qza_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/unweighted_unifrac_pcoa_results.qza'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'
        dokdo.beta_2d_plot(qza_file,
                           figsize=(5, 5))
        plt.tight_layout()

    .. code-block:: text

        # Explained proportions computed by QIIME 2:
        # 33.94% for Axis 1
        # 25.90% for Axis 2

    .. image:: images/beta_2d_plot-1.png

    We can color the datapoints with ``hue``. We can also change the
    style of datapoints with ``style``. If the variable of interest is
    numeric, we can use ``size`` to control the size of datapoints.
    Finally, we can combine all those groupings.

    .. code:: python3

        fig, [[ax1, ax2], [ax3, ax4]] = plt.subplots(2, 2, figsize=(15, 15))
        dokdo.beta_2d_plot(qza_file,
                           metadata_file,
                           ax=ax1,
                           hue='body-site')
        dokdo.beta_2d_plot(qza_file,
                           metadata_file,
                           ax=ax2,
                           style='subject')
        dokdo.beta_2d_plot(qza_file,
                           metadata_file,
                           ax=ax3,
                           size='days-since-experiment-start')
        dokdo.beta_2d_plot(qza_file,
                           metadata_file,
                           ax=ax4,
                           hue='body-site',
                           style='subject',
                           size='days-since-experiment-start')
        ax1.set_title("hue='body-site'", fontsize=20)
        ax2.set_title("style='subject'", fontsize=20)
        ax3.set_title("size='days-since-experiment-start'", fontsize=20)
        ax4.set_title('Multiple groupings', fontsize=20)
        for ax in [ax1, ax2, ax3, ax4]:
            ax.xaxis.label.set_size(20)
            ax.yaxis.label.set_size(20)
            ax.tick_params(axis='both', which='major', labelsize=15)
            ax.legend(loc='upper left')
        plt.tight_layout()

    .. image:: images/beta_2d_plot-2.png

    We can control categorical mapping of the ``hue`` variable with
    ``palette``:

    .. code:: python3

        palette = {'gut': 'yellow', 'left palm': 'green', 'right palm': 'blue', 'tongue': 'red'}
        dokdo.beta_2d_plot(qza_file, metadata_file, hue='body-site', palette=palette)
        plt.tight_layout()

    .. image:: images/beta_2d_plot-3.png
    """
    if isinstance(artifact, pd.DataFrame):
        df = artifact
    else:
        if isinstance(artifact, str):
            _pcoa_results = Artifact.load(artifact)
        else:
            _pcoa_results = artifact
        ordination_results = _pcoa_results.view(OrdinationResults)
        df = ordination_results.samples.iloc[:, :2]
        props = ordination_results.proportion_explained
        print('# Explained proportions computed by QIIME 2:')
        print(f'# {props[0]*100:.2f}% for Axis 1')
        print(f'# {props[1]*100:.2f}% for Axis 2')

    df.columns = ['Axis 1', 'Axis 2']

    if metadata is None:
        pass
    else:
        mf = common.get_mf(metadata)
        df = pd.concat([df, mf], axis=1, join='inner')

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    sns.scatterplot(
        x='Axis 1', y='Axis 2', data=df, hue=hue, hue_order=hue_order,
        style=style, style_order=style_order, size=size, ax=ax,
        s=s, alpha=alpha, legend=legend, palette=palette, **kwargs
    )

    return ax
