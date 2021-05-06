Dokdo API
*********

Introduction
============

This page describes Dokdo API, which is designed to be used with Jupyter Notebook in Python. Before using Dokdo API, please make sure your notebook is open within an environment where QIIME 2 and Dokdo are already installed.

First, at the beginning of your notebook, enter the following to import Dokdo API.

.. plot::
   :context:

   >>> import dokdo

Next, add the following to make figures. You should have the ``matplotlib`` package already installed in your environment because it is included in QIIME 2 installation. With the magic function ``%matplotlib inline``, the output of plotting methods will be displayed inline within Jupyter Notebook.

.. code:: python3

    import matplotlib.pyplot as plt
    %matplotlib inline

Finally, set the seed so that our results are reproducible.

.. code:: python3

    import numpy as np
    np.random.seed(1)

Tips
====

Setting Figure Properties
-------------------------

In this section, you'll learn how to control various properties of a figure using the plotting method ``dokdo.denoising_stats_plot`` as an example. This method creates a grouped box chart using denoising statistics from the DADA 2 algorithm. For more details about the method, see the ``dokdo.denoising_stats_plot`` section.

Let's start with a toy example. The figure below does not have a legend, which is bad, but let's not worry about that now.

.. plot::
    :context: close-figs

    >>> import matplotlib.pyplot as plt
    >>> data_dir = '/Users/sbslee/Desktop/dokdo/data/atacama-soil-microbiome-tutorial'
    >>> qza_file = f'{data_dir}/denoising-stats.qza'
    >>> metadata_file = f'{data_dir}/sample-metadata.tsv'
    >>> where = 'transect-name'
    >>> args = [qza_file, metadata_file, where]
    >>> dokdo.denoising_stats_plot(*args)
    >>> plt.tight_layout()

General Methods
===============

get_mf
------

.. automodule:: dokdo.api.get_mf
   :members:

ordinate
--------

.. automodule:: dokdo.api.ordinate
   :members:

pname
-----

.. automodule:: dokdo.api.pname
   :members:

num2sig
-------

.. automodule:: dokdo.api.num2sig
   :members:

wilcoxon
--------

.. automodule:: dokdo.api.wilcoxon
   :members:

mannwhitneyu
------------

.. automodule:: dokdo.api.mannwhitneyu
   :members:

Main Plotting Methods
=====================

read_quality_plot
-----------------

.. automodule:: dokdo.api.read_quality_plot
   :members:

denoising_stats_plot
--------------------

.. automodule:: dokdo.api.denoising_stats_plot
   :members:

alpha_rarefaction_plot
----------------------

.. automodule:: dokdo.api.alpha_rarefaction_plot
   :members:

alpha_diversity_plot
--------------------

.. automodule:: dokdo.api.alpha_diversity_plot
   :members:

beta_2d_plot
------------

.. automodule:: dokdo.api.beta_2d_plot
   :members:

beta_3d_plot
------------

.. automodule:: dokdo.api.beta_3d_plot
   :members:

beta_scree_plot
---------------

.. automodule:: dokdo.api.beta_scree_plot
   :members:

beta_parallel_plot
------------------

.. automodule:: dokdo.api.beta_parallel_plot
   :members:

distance_matrix_plot
--------------------

.. automodule:: dokdo.api.distance_matrix_plot
   :members:

taxa_abundance_bar_plot
-----------------------

.. automodule:: dokdo.api.taxa_abundance_bar_plot
   :members:

taxa_abundance_box_plot
-----------------------

.. automodule:: dokdo.api.taxa_abundance_box_plot
  :members:

ancom_volcano_plot
------------------

.. automodule:: dokdo.api.ancom_volcano_plot
  :members:

Other Plotting Methods
======================
