Dokdo CLI
*********

Introduction
============

This section describes command line interface (CLI) for the Dokdo package.

For getting help on CLI:

.. code-block:: console

    $ dokdo -h
    usage: dokdo [-v] [-h] COMMAND ...

    positional arguments:
      COMMAND        Name of the command.
        collapse     Create seven collapsed feature tables, one for each taxonomic
                     level (i.e. 'level-1.csv' to 'level-7.csv').
        make-manifest
                     Create a manifest file (.tsv) from a directory containing
                     FASTQ files.
        add-metadata
                     Add new metadata columns to an existing metadata file (.tsv).
        summarize    Extract summary or verbose data from an Artifact file.
        prepare-lefse
                     Create a TSV file which can be used as input for the LEfSe
                     tool.
        count-reads  Count the number of sequence reads from FASTQ.

    optional arguments:
      -v, --version  Show the version and exit.
      -h, --help     Show this help message and exit.

Commands
========

collapse
--------

.. code-block:: console

    $ dokdo collapse -h
    usage: dokdo collapse -t PATH -x PATH -o PATH [-h]

    Create seven collapsed feature tables, one for each taxonomic level (i.e.
    'level-1.csv' to 'level-7.csv').

    Arguments:
      -t PATH, --table-file PATH
                            Path to the table file with the
                            'FeatureTable[Frequency]' type. [required]
      -x PATH, --taxonomy-file PATH
                            Path to the taxonomy file with the
                            'FeatureData[Taxonomy]' type. [required]
      -o PATH, --output-dir PATH
                            Path to the output directory. [required]
      -h, --help            Show this help message and exit.

make-manifest
-------------

.. code-block:: console

    $ dokdo make-manifest -h
    usage: dokdo make-manifest -i PATH -o PATH [-h]

    Create a manifest file (.tsv) from a directory containing FASTQ files. The
    file names must include either '_R1_001.fastq' or '_R1_002.fastq'. The word
    before the third-to-last underscore will be set as the sample ID. For example,
    a file named 'EXAMPLE_S1_R1_001.fastq.gz' will produce 'EXAMPLE' as sample ID
    and 'EXAM_PLE_S1_R1_001.fastq.gz', 'EXAM_PLE'.

    Arguments:
      -i PATH, --fastq-dir PATH
                            Path to the directory containing input FASTQ files.
                            [required]
      -o PATH, --output-file PATH
                            Path to the output file. [required]
      -h, --help            Show this help message and exit.

add-metadata
------------

.. code-block:: console

    $ dokdo add-metadata -h
    usage: dokdo add-metadata [-i PATH] [-c PATH] [-o PATH] [-h]

    Add new metadata columns to an existing metadata file (.tsv). The files
    '-i/--metadata-file' and '-c/--columns-file' must have at least one
    overlapping column.

    Arguments:
      -i PATH, --metadata-file PATH
                            Path to the metadata file. [required]
      -c PATH, --columns-file PATH
                            Path to a text file (.tsv) containing the columns to
                            be added. The first row should be column names.
                            [required]
      -o PATH, --output-file PATH
                            Path to the output file. [required]
      -h, --help            Show this help message and exit.

summarize
---------

.. code-block:: console

    $ dokdo summarize -h
    usage: dokdo summarize [-i PATH] [-v] [-h]

    Extract summary or verbose data from an Artifact file. This command
    automatically detects the input file's semantic type and then extracts summary
    or verbose data from it. Currently, the command supports the following
    semantic types: FeatureTable[Frequency], FeatureTable[RelativeFrequency],
    FeatureData[Sequence], FeatureData[AlignedSequence], FeatureData[Taxonomy],
    DistanceMatrix.

    Arguments:
      -i PATH, --input-file PATH
                            Path to the input Artifact file. [required]
      -v, --verbose         Print a verbose version of the results.
      -h, --help            Show this help message and exit.

prepare-lefse
-------------

.. code-block:: console

    $ dokdo prepare-lefse -h
    usage: dokdo prepare-lefse -t PATH -x PATH -m PATH -o PATH -c TEXT [-s TEXT]
                               [-u TEXT] [-w TEXT] [-h]

    Create a TSV file which can be used as input for the LEfSe tool. This command
    1) collapses the input feature table at the genus level, 2) computes relative
    frequency of the features, 3) performs sample filtration if requested, 4)
    changes the format of feature names, 5) adds the relevant metadata as 'Class',
    'Subclass', and 'Subject', and 6) writes a text file which can be used as
    input for LEfSe.

    Arguments:
      -t PATH, --table-file PATH
                            Path to the table file with the
                            'FeatureTable[Frequency]' type. [required]
      -x PATH, --taxonomy-file PATH
                            Path to the taxonomy file with the
                            'FeatureData[Taxonomy]' type. [required]
      -m PATH, --metadata-file PATH
                            Path to the metadata file. [required]
      -o PATH, --output-file PATH
                            Path to the output file. [required]
      -c TEXT, --class-col TEXT
                            Metadata column used as 'Class' by LEfSe. [required]
      -s TEXT, --subclass-col TEXT
                            Metadata column used as 'Subclass' by LEfSe.
      -u TEXT, --subject-col TEXT
                            Metadata column used as 'Subject' by LEfSe.
      -w TEXT, --where TEXT
                            SQLite 'WHERE' clause specifying sample metadata
                            criteria.
      -h, --help            Show this help message and exit.
