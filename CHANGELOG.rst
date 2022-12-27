Changelog
*********

1.16.0 (2022-12-27)
-------------------

* Add new argument ``palette`` to :meth:`taxa_abundance_box_plot` method.
* :issue:`51`: Fix bug in :meth:`clustermap` method where color labels disappeared when ``flip=True``.

1.15.0 (2022-08-23)
-------------------

* Update :meth:`pname` method to support taxa names that don't have a delimiter such as ASV IDs (e.g. ``1ad289cd8f44e109fd95de0382c5b252``). Basically, the method will simply return the input name as is.
* Add new argument ``delimiter`` to :meth:`pname` method.
* Add new method :meth:`group_correlation_heatmap`.

1.14.0 (2022-06-24)
-------------------

* :issue:`39`: Add new argument ``palette`` to :meth:`beta_2d_plot` and :meth:`beta_3d_plot` methods so that users can control categorical color mapping.
* :issue:`39`: Bring back ``add_datapoints`` argument to :meth:`taxa_abundance_box_plot` method to display data points on top of the boxes (note that it was deprecated in 1.11.0 version).
* Add new argument ``kwargs`` to :meth:`beta_2d_plot` method. It will be passed to :meth:`seaborn.scatterplot()` method so that users can have more control over various aspects of the output plot.
* :issue:`41`: Fix bug in :meth:`beta_parallel_plot` method when users provide the ``ax`` argument.

1.13.0 (2022-04-02)
-------------------

* :issue:`37`: Add new methods :meth:`cross_association_table`, :meth:`cross_association_heatmap`, and :meth:`cross_association_regplot`.

1.12.0 (2022-02-11)
-------------------

* :issue:`33`: Update the :command:`make-manifest` command to ignore undetermined FASTQ files.
* :issue:`34`: Update the :meth:`alpha_rarefaction_plot` method to keep 'N/A' values as string instead of NaN.
* :issue:`35`: Update the methods :meth:`alpha_diversity_plot`, :meth:`beta_2d_plot`, and :meth:`beta_3d_plot` to accept :class:`pandas.DataFrame` in case the input data was not generated from QIIME 2 (e.g. shotgun sequencing).
* Update the methods :meth:`beta_2d_plot` and :meth:`beta_3d_plot` to print out the proportions explained instead of embedding them in the PCoA plot.

1.11.0 (2021-08-04)
-------------------

* :issue:`21`: Update the :meth:`taxa_abundance_bar_plot` method to accept :class:`pandas.DataFrame` in case the input data was not generated from QIIME 2 (e.g. shotgun sequencing).
* :issue:`19`: Remove the ``artist_kwargs`` argument from all the remaining methods:

  - :meth:`denoising_stats_plot`
  - :meth:`alpha_rarefaction_plot`
  - :meth:`alpha_diversity_plot`
  - :meth:`beta_2d_plot`
  - :meth:`beta_3d_plot`
  - :meth:`beta_scree_plot`
  - :meth:`beta_parallel_plot`
  - :meth:`distance_matrix_plot`
  - :meth:`regplot`
  - :meth:`taxa_abundance_box_plot`
  - :meth:`taxa_abundance_bar_plot`

* Deprecate the :meth:`barplot` method.
* :issue:`22`: Rename the :meth:`heatmap` method to :meth:`clustermap`.
* :issue:`22`: Update the :meth:`clustermap` method to accept :class:`pandas.DataFrame` in case the input data was not generated from QIIME 2 (e.g. shotgun sequencing). You can now also flip the x and y axes with the ``flip`` option.
* :issue:`22`: Add a new main plotting method :meth:`heatmap`.
* :issue:`24`: Update the :meth:`pname` method to allow returning of more than one tax level.
* Deprecate the ``add_datapoints`` argument in the :meth:`taxa_abundance_box_plot` method.

1.10.0 (2021-07-06)
-------------------

* :issue:`14`, :issue:`17`: Add the ``group_order`` option to the :meth:`taxa_abundance_bar_plot` method.
* Fix a minor bug in the :meth:`addbiplot` method when feature is 'Unassigned'.
* Deprecate the :command:`count-reads` command.
* :issue:`19`: Remove the ``artist_kwargs`` argument from the following methods:

  - :meth:`ancom_volcano_plot`
  - :meth:`read_quality_plot`

1.9.0 (2021-06-07)
------------------

* Add publicly available datasets from QIIME 2 for tutorials.
* :issue:`14`: Add the ``group`` option to the :meth:`taxa_abundance_bar_plot` method. Using this option will create a bar for each group instead of each sample.

1.8.0 (2021-05-09)
------------------

* Updated docstring.
* Moved the official documentation from Wiki page to Read the Docs.

1.7.0 (2021-04-05)
------------------

- Added a new command called ``count-reads`` which counts the number of sequence reads from FASTQ.
- Updated the ``summarize`` command.
- Updated the following methods:

    - ``taxa_abundance_box_plot()``
    - ``taxa_abundance_bar_plot()``
    - ``distance_matrix_plot()``
    - ``ordinate()``
    - ``barplot()``

- See :issue:`10` for more details.

1.6.0 (2021-03-08)
------------------

- Added a new method called ``pname()`` which returns a prettified taxon name.
- Added a new method called ``num2sig()`` which converts a p-value to significance annotation.
- Added a new method called ``wilcoxon()`` which performs the Wilcoxon Signed-rank test between two paired groups for a given taxon.
- Added a new method called ``mannwhitneyu()`` which performs the Mann–Whitney U test between two groups for a given taxon.
- There have been major changes to the ``heatmap()`` method. First, it now supports two grouping variables instead of just one (e.g. ``hue1`` and ``hue2``). Second, it supports the centered log-ratio (CLR) transformation as a normalization option (in addition to ``log10``). Third, it now has ``kwargs`` that are passed to the ``seaborn.clustermap()`` method (e.g. ``xticklabels=False``). Fourth, the bug giving the ``FloatingPointError: NaN dissimilarity value.`` error when sample-filtered metadata is provided and the ``metric='correlation'`` argument is used has been fixed. Fifth, the bug giving an error when one of the metadata columns has only zeros has been fixed.
- In addition to ``heatmap()``, the following methods have been updated:

    - ``addpairs()``
    - ``alpha_diversity_plot()``

- Updated the ``summarize`` command.
- Updated the ``prepare-lefse`` command to output more informative taxa name than just underscores (e.g. ``__`` and ``g__``).
- See :issue:`8` for more details.

1.5.0 (2021-02-03)
------------------

- Starting this version, Dokdo is packaged with ``setuptools``.
- There have been major changes to Dokdo CLI.
- Added a new plotting method called ``regplot()``.
- Added a new command called ``prepare-lefse``.
- The ``merge_metadata`` command has been deprecated.
- Updated the following methods:

    - ``_artist()``
    - ``alpha_diversity_plot()``
    - ``beta_3d_plot()``
    - ``beta_parallel_plot()``
    - ``barplot()``
    - ``ordinate()``
    - ``taxa_abundance_bar_plot()``
    - ``taxa_abundance_box_plot()``
    - ``heatmap()``

- Updated the ``make_manifest`` command.
- See :issue:`6` for more details.

1.4.0 (2021-01-09)
------------------

- Added a new command called ``summarize``.
- Added a new plotting method called ``heatmap()``.
- Updated the following commands:

    - ``make_manifest``
    - ``add_metadata``
    - ``collapse``

* Updated the following methods:

    - ``_artist()``
    - ``alpha_rarefaction_plot()``
    - ``taxa_abundance_bar_plot()``
    - ``taxa_abundance_box_plot()``

- See :issue:`4` for more details.

1.3.0 (2020-12-23)
------------------

- Updated the ``ordinate()`` method so that the user can now choose to:

    - skip rarefying,
    - provide custom sampling depth for rarefying,
    - provide ``qiime2.Artifact`` as input instead of file path, and
    - output ``PCoAResults % Properties('biplot')`` as well as ``PCoAResults``.

- Added new plotting methods:

    - ``beta_scree_plot()``
    - ``beta_parallel_plot()``
    - ``addbiplot()``
    - ``barplot()``

- See :issue:`2` for more details.

1.2.0 (2020-12-08)
------------------

- The ``tax2seq`` command has been deprecated.
- Updated the ``_artist()`` method to set the font size of title, labels, etc.
- Added the ``s`` argument to the ``ancom_volcano_plot()`` method for setting marker size.
- Updated the docstring.
- See :issue:`1` for more details.

1.1.0 (2020-11-23)
------------------

- Introduced the ``addpairs()`` method.
- The ``beta_2d_plot_gallery()`` method has been deprecated.
- Made some changes to the following methods:

    - ``ordinate()``
    - ``taxa_abundance_bar_plot()``
    - ``taxa_abundance_box_plot()``
    - ``_artist()``

- Fixed some bugs.
- Made keyword arguments for the ``_artist()`` method more explicit with ``artist_kwargs``.
- Temporary files will be deleted automatically from now on.
- Updated the docstring.
- Plotting methods now accept Artifact and Visualization objects as input.

1.0.0 (2020-11-09)
------------------

- Initial release.
