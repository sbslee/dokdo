QIIME 2 API
***********

Introduction
============

This page describes QIIME 2 API, which is the main building block of Dokdo, and how I use it to perform various tasks in Python. If you are planning to use QIIME 2 API with Jupyter Notebook, make sure your notebook is open within an environment with QIIME 2.

Get Sample IDs from Metadata
============================

.. code:: python3

    >>> from qiime2 import Metadata
    >>> metadata = Metadata.load('/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv')
    >>> ids = metadata.get_ids("[body-site] IN ('gut', 'right palm') AND [subject]='subject-1'")
    >>> print(list(ids)) # Conver set to list.
    ['L3S242', 'L1S57', 'L3S360', 'L3S294', 'L1S76', 'L3S341', 'L1S105', 'L3S313', 'L1S8']

Filter Samples from ASV Table
=============================

.. code:: python3

    >>> from qiime2 import Artifact
    >>> from qiime2 import Metadata
    >>> from qiime2.plugins import feature_table
    >>> table = Artifact.load('/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/table.qza')
    >>> metadata = Metadata.load('/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv')
    >>> filter_result = feature_table.methods.filter_samples(table=table, metadata=metadata, where="[subject]='subject-1'")
    >>> filtered_table = filter_result.filtered_table

Rarefy ASV Table
================

.. code:: python3

    >>> from qiime2 import Artifact
    >>> from qiime2.plugins import feature_table
    >>> table = Artifact.load('/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/table.qza')
    >>> rarefy_result = feature_table.methods.rarefy(table=table, sampling_depth=1103)
    >>> rarefy_result.rarefied_table.save('rarefied-table.qza')

Compute Distance Matrix
=======================

.. code:: python3

    >>> from qiime2 import Artifact
    >>> from qiime2.plugins import diversity_lib
    >>> rarefied_table = Artifact.load('output/QIIME-2-API/rarefied-table.qza')
    >>> phylogeny = Artifact.load('data/moving-pictures-tutorial/rooted-tree.qza')
    >>> distance_matrix_result = diversity_lib.methods.jaccard(table=rarefied_table)
    >>> # distance_matrix_result = diversity_lib.methods.bray_curtis(table=rarefied_table)
    >>> # distance_matrix_result = diversity_lib.methods.unweighted_unifrac(table=rarefied_table, phylogeny=phylogeny)
    >>> # distance_matrix_result = diversity_lib.methods.weighted_unifrac(table=rarefied_table, phylogeny=phylogeny)
    >>> distance_matrix = distance_matrix_result.distance_matrix
    >>> distance_matrix.save('jaccard_distance_matrix.qza')

Visualize Ordination
====================

.. code:: python3

    >>> from qiime2 import Artifact
    >>> from qiime2 import Metadata
    >>> from qiime2.plugins import diversity
    >>> from qiime2.plugins import emperor
    >>> distance_matrix = Artifact.load('data/moving-pictures-tutorial/jaccard_distance_matrix.qza')
    >>> metadata = Metadata.load('data/moving-pictures-tutorial/sample-metadata.tsv')
    >>> pcoa_result = diversity.methods.pcoa(distance_matrix=distance_matrix)
    >>> pcoa = pcoa_result.pcoa
    >>> pcoa.save('output/QIIME-2-API/jaccard_pcoa_results.qza')
    >>> emperor_result = emperor.visualizers.plot(pcoa=pcoa, metadata=metadata)
    >>> emperor = emperor_result.visualization
    >>> emperor.save('output/QIIME-2-API/jaccard_emperor.qzv')
    >>> emperor

Alpha Rarefaction
=================

.. code:: python3

    from qiime2 import Artifact
    from qiime2 import Metadata
    from qiime2.plugins import diversity
    >>> alpha_rarefaction_result = diversity.visualizers.alpha_rarefaction(table=Artifact.load('data/moving-pictures-tutorial/table.qza'),
    ...                                                                    phylogeny=Artifact.load('data/moving-pictures-tutorial/rooted-tree.qza'),
    ...                                                                    max_depth=4000,
    ...                                                                    metadata=Metadata.load('data/moving-pictures-tutorial/sample-metadata.tsv'))
    >>> rarefaction = alpha_rarefaction_result.visualization
    >>> rarefaction.save('output/QIIME-2-API/alpha-rarefaction.qzv')
    >>> rarefaction

Note that if you do not include the ``metadata`` option, the visualization will only show sample-level results.

Collapse Feature Table
======================

.. code:: python3

    >>> from qiime2 import Artifact
    >>> from qiime2.plugins import taxa
    >>> table_file = 'data/moving-pictures-tutorial/table.qza'
    >>> taxonomy_file = 'data/moving-pictures-tutorial/taxonomy.qza'
    >>> collapse_result = taxa.methods.collapse(table=Artifact.load(table_file),
    ...                                         taxonomy=Artifact.load(taxonomy_file),
    ...                                         level=6)
    >>> collapsed_table = collapse_result.collapsed_table
