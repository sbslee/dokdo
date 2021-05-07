Changelog
*********

1.8.0 (in development)
----------------------

* Transitioned to Read the Docs

1.7.0 (2021-04-05)
------------------

* Added a new command called `count-reads` which counts the number of sequence reads from FASTQ.
* Updated the `summarize` command.
* Updated the following methods:
    * `taxa_abundance_box_plot()`
    * `taxa_abundance_bar_plot()`
    * `distance_matrix_plot()`
    * `ordinate()`
    * `barplot()`
* See :issue:`10` for more details.

1.6.0 (2021-03-08)
------------------

* Added a new method called `pname()` which returns a prettified taxon name.
* Added a new method called `num2sig()` which converts a p-value to significance annotation.
* Added a new method called `wilcoxon()` which performs the Wilcoxon Signed-rank test between two paired groups for a given taxon.
* Added a new method called `mannwhitneyu()` which performs the Mann–Whitney U test between two groups for a given taxon.
* There have been major changes to the `heatmap()` method. First, it now supports two grouping variables instead of just one (e.g. `hue1` and `hue2`). Second, it supports the centered log-ratio (CLR) transformation as a normalization option (in addition to `log10`). Third, it now has `kwargs` that are passed to the `seaborn.clustermap()` method (e.g. `xticklabels=False`). Fourth, the bug giving the `FloatingPointError: NaN dissimilarity value.` error when sample-filtered metadata is provided and the `metric='correlation'` argument is used has been fixed. Fifth, the bug giving an error when one of the metadata columns has only zeros has been fixed.
* In addition to `heatmap()`, the following methods have been updated:
    * `addpairs()`
    * `alpha_diversity_plot()`
* Updated the `summarize` command.
* Updated the `prepare-lefse` command to output more informative taxa name than just underscores (e.g. `__` and `g__`).
* See :issue:`8` for more details.

1.5.0 (2021-02-03)
------------------

* Starting this version, Dokdo is packaged with `setuptools`.
* There have been major changes to the Dokdo CLI.
* Added a new plotting method called `regplot()`.
* Added a new command called `prepare-lefse`.
* The `merge_metadata` command has been deprecated.
* Updated the following methods:
    * `_artist()`
    * `alpha_diversity_plot()`
    * `beta_3d_plot()`
    * `beta_parallel_plot()`
    * `barplot()`
    * `ordinate()`
    * `taxa_abundance_bar_plot()`
    * `taxa_abundance_box_plot()`
    * `heatmap()`
* Updated the `make_manifest` command.
* See :issue:`6` for more details.

1.4.0 (2021-01-09)
------------------

* Added a new command called `summarize`.
* Added a new plotting method called `heatmap()`.
* Updated the following commands:
    * `make_manifest`
    * `add_metadata`
    * `collapse`
* Updated the following methods:
    * `_artist()`
    * `alpha_rarefaction_plot()`
    * `taxa_abundance_bar_plot()`
    * `taxa_abundance_box_plot()`
* See :issue:`4` for more details.

1.3.0 (2020-12-23)
------------------

* Updated the `ordinate()` method so that the user can now choose to:
    * skip rarefying,
    * provide custom sampling depth for rarefying,
    * provide `qiime2.Artifact` as input instead of file path, and
    * output `PCoAResults % Properties('biplot')` as well as `PCoAResults`.
* Added new plotting methods:
    * `beta_scree_plot()`
    * `beta_parallel_plot()`
    * `addbiplot()`
    * `barplot()`
* See :issue:`2` for more details.

1.2.0 (2020-12-08)
------------------

* The `tax2seq` command has been deprecated.
* Updated the `_artist()` method to set the font size of title, labels, etc.
* Added the `s` argument to the `ancom_volcano_plot()` method for setting marker size.
* Updated the docstring.
* See :issue:`1` for more details.

1.1.0 (2020-11-23)
------------------

* Introduced the `addpairs()` method.
* The `beta_2d_plot_gallery()` method has been deprecated.
* Made some changes to the following methods:
    * `ordinate()`
    * `taxa_abundance_bar_plot()`
    * `taxa_abundance_box_plot()`
    * `_artist()`
* Fixed some bugs.
* Made keyword arguments for the `_artist()` method more explicit with `artist_kwargs`.
* Temporary files will be deleted automatically from now on.
* Updated the docstring.
* Plotting methods now accept Artifact and Visualization objects as input.

1.0.0 (2020-11-09)
------------------

* Initial release.
