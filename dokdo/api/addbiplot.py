from skbio.stats.ordination import OrdinationResults
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean
from qiime2 import Artifact

def addbiplot(
    pcoa_results, dim=2, scale=1.0, count=5, fontsize=None,
    name_type='feature', taxonomy=None, level=None, ax=None, figsize=None
):
    """
    Draw arrows (features) to an existing PCoA plot.

    This method supports both 2D and 3D plots.

    Parameters
    ----------
    pcoa_results : str or qiime2.Artifact
        Artifact file or object corresponding to
        PCoAResults % Properties('biplot').
    dim : [2, 3], default: 2
        Dimension of the plot.
    scale : float, default: 1.0
        Scale for arrow length.
    count : int, default: 5
        Number of important features to be displayed.
    fontsize : float or str, optional
        Sets font size.
    name_type : ['feature', 'taxon', 'confidence'], default: 'feature'
        Determines the type of names displayed. Using 'taxon' and 'confidence'
        requires taxonomy.
    taxonomy : str or qiime2.Artifact
        Artifact file or object corresponding to FeatureData[Taxonomy].
        Required if ``name_type`` is 'taxon' or 'confidence'.
    level : int, optional
        Level of taxonomic rank to be displayed. This argument has an effect
        only when ``name_type`` is 'taxon'.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    dokdo.api.ordinate
    dokdo.api.beta_2d_plot
    dokdo.api.beta_3d_plot

    Examples
    --------
    Below is a simple example.

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()

        table_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/table.qza'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'
        pcoa_results = dokdo.ordinate(table_file, sampling_depth=0, biplot=True, number_of_dimensions=10)
        ax = dokdo.beta_2d_plot(pcoa_results,
                                hue='body-site',
                                metadata=metadata_file,
                                figsize=(8, 8))
        dokdo.addbiplot(pcoa_results, ax=ax, count=3)
        plt.tight_layout()

    .. image:: images/addbiplot-1.png

    We can also draw a 3D biplot.

    .. code:: python3

        ax = dokdo.beta_3d_plot(pcoa_results,
                                hue='body-site',
                                metadata=metadata_file,
                                figsize=(8, 8))
        dokdo.addbiplot(pcoa_results, ax=ax, count=3, dim=3)
        plt.tight_layout()

    .. image:: images/addbiplot-2.png

    Finally, we can display taxonomic classification instead of feature ID.

    .. code:: python3

        taxonomy_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/taxonomy.qza'
        ax = dokdo.beta_3d_plot(pcoa_results,
                                hue='body-site',
                                metadata=metadata_file,
                                figsize=(8, 8))
        dokdo.addbiplot(pcoa_results,
                        ax=ax,
                        count=3,
                        dim=3,
                        taxonomy=taxonomy_file,
                        name_type='taxon',
                        level=6)
        plt.tight_layout()

    .. image:: images/addbiplot-3.png
    """
    if isinstance(pcoa_results, str):
        _pcoa_results = Artifact.load(pcoa_results)
    else:
        _pcoa_results = pcoa_results

    ordination_results = _pcoa_results.view(OrdinationResults)

    feats = ordination_results.features.copy()
    origin = np.zeros_like(feats.columns)
    feats['importance'] = feats.apply(euclidean, axis=1, args=(origin,))
    feats.sort_values('importance', inplace=True, ascending=False)
    feats.drop(['importance'], inplace=True, axis=1)
    feats = feats[:count]

    if taxonomy is not None:
        if isinstance(taxonomy, str):
            _taxonomy = Artifact.load(taxonomy)
        else:
            _taxonomy = taxonomy
        tax_df = _taxonomy.view(pd.DataFrame)
        feats = pd.concat([feats, tax_df], axis=1, join='inner')

    def f(s):
        if s == 'Unassigned':
            return s
        ranks = list(s.split(';'))
        if level is None:
            return s
        else:
            return ranks[level-1]

    if name_type == 'feature':
        names = feats.index
    elif name_type == 'taxon':
        names = [f(x) for x in feats['Taxon']]
    else:
        names = feats['Confidence']

    if ax is None:
        if dim == 2:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = plt.figure(figsize=figsize)
            ax = fig.add_subplot(1, 1, 1, projection='3d')

    for i in range(len(feats)):
        x = feats.iloc[i, 0]*scale
        y = feats.iloc[i, 1]*scale
        z = feats.iloc[i, 2]*scale

        a = [[0, x], [0, y]]
        b = [x, y]

        if dim == 3:
            a.append([0, z])
            b.append(z)

        ax.plot(*a, color='black')
        ax.text(*b, names[i], ha='center', fontsize=fontsize)

    return ax
