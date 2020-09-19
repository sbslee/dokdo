# dokdo

```dokdo``` is a collection of scripts for microbiome sequencing analysis. For details, please see the Wiki page. Pull requests are welcome.

# Command Line Interface (CLI)

```
$ python3 /path/to/dokdo/cli.py -h
usage: cli.py [-h] command ...

positional arguments:
  command     name of command
    collapse  create 7 collapsed ASV tables, one for each taxonomic level
    tax2seq   return mapping between ASVs and taxonomic classifications

optional arguments:
  -h, --help  show this help message and exit
```