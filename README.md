# dokdo

```dokdo``` is a collection of scripts for microbiome sequencing analysis. For details, please see the Wiki page. Pull requests are welcome.

# Command Line Interface (CLI)

```
$ python3 /path/to/dokdo/cli.py -h
usage: cli.py [-h] command ...

positional arguments:
  command        Name of the command.
    collapse     This command creates seven collapsed ASV tables, one for each
                 taxonomic level.
    tax2seq      This command returns the mapping between observed ASVs and
                 taxonomic classifications.
    make_manifest
                 This command creates a manifest file from a directory
                 containing FASTQ files.
    add_metadata
                 This command adds new columns to an existing sample-metadata
                 file.

optional arguments:
  -h, --help     show this help message and exit
```