README
******

.. image:: https://readthedocs.org/projects/dokdo/badge/?version=latest
   :target: https://dokdo.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

Dokdo is a lightweight Python package for microbiome sequencing analysis, which can be used as a command line tool and as a Python module. Dokdo is designed to be used with `QIIME 2 <https://qiime2.org/>`_, a powerful community-developed platform for microbiome bioinformatics.

Dokdo is versatile like the swiss army knife. For example, you can use its command line interface to create various input files for QIIME 2 (e.g. a manifest file or a sample-metadata file) or perform a variety of secondary analyses. You may also use Dokdo's application programming interface with `Jupyter Notebook <https://jupyter.org/>`_ to create publication-quality figures using the output files from QIIME 2 (e.g. a taxonomic bar plot or an alpha rarefaction plot).

To install Dokdo, enter the following in your terminal:

.. code-block:: console

   $ git clone https://github.com/sbslee/dokdo
   $ cd dokdo
   $ pip install .

For more details, please see the `Read the Docs <https://dokdo.readthedocs.io/en/latest/>`_ which serves as the official user documentation for Dokdo, including instructions, tutorials, and other important information. The page also documents select information and resources pertaining to QIIME 2 that I find particularly useful.

Pull requests are always welcome.

| Author: Seung-been "Steven" Lee
| Email: sbstevenlee@gmail.com
| License: MIT License
