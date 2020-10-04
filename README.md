# dokdo

## Table of Contents

* [Introduction](#Introduction)
* [Wiki](#Wiki)
* [Command Line Interface](#Command-Line-Interface)
* [Application Programming Interface](#Application-Programming-Interface)


## Introduction <a name="Introduction"></a>

Dokdo is a Python package for microbiome sequencing analysis, which can be used as a standalone program and as a Python module. Dokdo is designed to be used with [QIIME 2](https://qiime2.org/). For details, please see the Wiki page. Pull requests are welcome.

## Wiki <a name="Wiki"></a>

The Wiki page is the official user documentation for Dokdo, including instructions, tutorials, and other important information. The page also documents select information and resources pertaining to QIIME 2 that I find particularly useful.

## Command Line Interface <a name="Command-Line-Interface"></a>

For getting help:

```
$ python3 /path/to/dokdo/cli.py -h
usage: cli.py [-h] command ...

positional arguments:
  command         Name of the command.
    collapse      This command creates seven collapsed ASV tables, one for
                  each taxonomic level.
    tax2seq       This command returns the mapping between observed ASVs and
                  taxonomic classifications.
    make_manifest
                  This command creates a manifest file from a directory
                  containing FASTQ files.
    add_metadata  This command adds new columns to an existing sample-metadata
                  file.
    merge_metadata
                  This command merges two or more sample-metadata.tsv files
                  vertically.

optional arguments:
  -h, --help      show this help message and exit
```

For getting command-specific help:

```
$ python3 /path/to/dokdo/cli.py add_metadata -h
usage: cli.py add_metadata [-h] metadata columns output

This command adds new columns to an existing sample-metadata file. The
'metadata' file and the 'columns' file must have at least one overlapping
column.

positional arguments:
  metadata    Path to the input sample-metadata file (.tsv).
  columns     Path to a file containing the new columns to be added (.tsv).
              The first row should be column names.
  output      Path to the output sample-metadata file (.tsv).

optional arguments:
  -h, --help  show this help message and exit
```

## Application Programming Interface <a name="Application-Programming-Interface"></a>

At the beginning of your Jupyter Notebook, enter the following:

```
import sys
sys.path.append("/path/to/dokdo")
import api
```

Then you can use the methods defined in Dokdo API. For example:

```
help(api.ancom_volcano_plot)
```

Gives:

```
Help on function ancom_volcano_plot in module api:

ancom_volcano_plot(ancom, ax=None)
    This method creates an ANCOM volcano plot.
    
    Example:
        import matplotlib.pyplot as plt
        api.ancom_volcano_plot('ancom-Site.qzv')
        plt.savefig('ancom-Site.pdf')
```