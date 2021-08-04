from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins import diversity_lib
from qiime2.plugins import feature_table
from qiime2.plugins import diversity
import pandas as pd

def ordinate(
    table, metadata=None, metric='jaccard', sampling_depth=-1,
    phylogeny=None, number_of_dimensions=None, biplot=False
):
    """
    Perform ordination using principal coordinate analysis (PCoA).

    This method wraps multiple QIIME 2 methods to perform ordination and
    returns Artifact object containing PCoA results.

    Under the hood, this method filters the samples (if requested), performs
    rarefying of the feature table (if requested), computes distance matrix,
    and then runs PCoA.

    By default, the method returns PCoAResults. For creating a biplot,
    use `biplot=True` which returns PCoAResults % Properties('biplot').

    Parameters
    ----------
    table : str or qiime2.Artifact
        Artifact file or object corresponding to FeatureTable[Frequency].
    metadata : str or qiime2.Metadata, optional
        Metadata file or object. All samples in 'metadata' that are also in
        the feature table will be retained.
    metric : str, default: 'jaccard'
        Metric used for distance matrix computation ('jaccard',
        'bray_curtis', 'unweighted_unifrac', or 'weighted_unifrac').
    sampling_depth : int, default: -1
        If negative, skip rarefying. If 0, rarefy to the sample with minimum
        depth. Otherwise, rarefy to the provided sampling depth.
    phylogeny : str, optional
        Rooted tree file. Required if using 'unweighted_unifrac', or
        'weighted_unifrac' as metric.
    number_of_dimensions : int, optional
        Dimensions to reduce the distance matrix to.
    biplot : bool, default: False
        If true, return PCoAResults % Properties('biplot').

    Returns
    -------
    qiime2.Artifact
        Artifact object corresponding to PCoAResults or
        PCoAResults % Properties('biplot').

    See Also
    --------
    dokdo.api.beta_2d_plot
    dokdo.api.beta_3d_plot
    dokdo.api.beta_scree_plot
    dokdo.api.beta_parallel_plot

    Notes
    -----
    The resulting Artifact object can be directly used for plotting.

    Examples
    --------
    Below is a simple example. Note that the default distance metric
    used is ``jaccard``. The resulting object ``pcoa`` can be directly
    used for plotting by the ``dokdo.beta_2d_plot`` method as shown below.

    .. code:: python3

        import dokdo
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()

        qza_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/table.qza'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'

        pcoa_results = dokdo.ordinate(qza_file)

        dokdo.beta_2d_plot(
            pcoa_results,
            metadata=metadata_file,
            hue='body-site',
            figsize=(8, 8)
        )

        plt.tight_layout()

    .. image:: images/ordinate-1.png

    You can choose a subset of samples.

    .. code:: python3

        from qiime2 import Metadata

        mf = dokdo.get_mf(metadata_file)
        mf = mf[mf['body-site'].isin(['gut', 'left palm'])]

        pcoa_results = dokdo.ordinate(qza_file, metadata=Metadata(mf))

        dokdo.beta_2d_plot(
            pcoa_results,
            metadata=metadata_file,
            hue='body-site',
            figsize=(8, 8)
        )

        plt.tight_layout()

    .. image:: images/ordinate-2.png

    You can also generate a biplot.

    .. code:: python3

        pcoa_results = dokdo.ordinate(qza_file, biplot=True, number_of_dimensions=10)

        ax = dokdo.beta_2d_plot(
            pcoa_results,
            metadata=metadata_file,
            hue='body-site',
            figsize=(8, 8)
        )


        dokdo.addbiplot(pcoa_results, ax=ax, count=7)

        plt.tight_layout()

    .. image:: images/ordinate-3.png
    """
    if isinstance(table, Artifact):
        table = table
    elif isinstance(table, str):
        table = Artifact.load(table)
    else:
        raise TypeError(f"Incorrect feature table type: {type(table)}")

    # If metadata is provided, perform sample filtration.
    if metadata is not None:
        if isinstance(metadata, Metadata):
            _metadata = metadata
        else:
            _metadata = Metadata.load(metadata)
        _table = feature_table.methods.filter_samples(
            table=table, metadata=_metadata).filtered_table
    else:
        _table = table

    # Perform rarefying.
    if sampling_depth < 0:
        rarefied_table = _table
    else:
        if sampling_depth == 0:
            sampling_depth = int(_table.view(pd.DataFrame).sum(axis=1).min())

        rarefy_result = feature_table.methods.rarefy(
            table=_table, sampling_depth=sampling_depth)

        rarefied_table = rarefy_result.rarefied_table

    if metric == 'jaccard':
        distance_matrix_result = diversity_lib.methods.jaccard(
            table=rarefied_table)
    elif metric == 'bray_curtis':
        distance_matrix_result = diversity_lib.methods.bray_curtis(
            table=rarefied_table)
    elif metric == 'unweighted_unifrac':
        distance_matrix_result = diversity_lib.methods.unweighted_unifrac(
            table=rarefied_table, phylogeny=Artifact.load(phylogeny))
    elif metric == 'weighted_unifrac':
        distance_matrix_result = diversity_lib.methods.weighted_unifrac(
            table=rarefied_table, phylogeny=Artifact.load(phylogeny))
    else:
        raise ValueError(f"Incorrect metric detected: {metric}")

    distance_matrix = distance_matrix_result.distance_matrix

    result_obj = diversity.methods.pcoa(distance_matrix=distance_matrix,
        number_of_dimensions=number_of_dimensions)
    pcoa_results = result_obj.pcoa

    if biplot:
        rf_result = feature_table.methods.relative_frequency(table=_table)
        rf_table = rf_result.relative_frequency_table
        result_obj = diversity.methods.pcoa_biplot(pcoa=pcoa_results,
            features=rf_table)
        pcoa_results = result_obj.biplot

    return pcoa_results
