import argparse

import pandas as pd

from qiime2 import Artifact

from qiime2.plugins import taxa

def collapse(table, taxonomy, output_dir):
    for i in range(1, 8):
        _ = taxa.methods.collapse(table=Artifact.load(table),
                                  taxonomy=Artifact.load(taxonomy),
                                  level=i)
        _.collapsed_table.view(pd.DataFrame).T.to_csv(f'level-{i}.csv')

def main():
    commands = {
        'collapse': collapse,
    }

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(dest='command',
                                       metavar='command',
                                       help="name of command")
    subparsers.required = True

    collapse_parser = subparsers.add_parser(
        'collapse',
        help="create 7 collapsed ASV tables, one for each taxonomic level",
    )
    collapse_parser.add_argument('table',
                                 help="path to input table.qza file")
    collapse_parser.add_argument('taxonomy',
                                 help="path to input taxonomy.qza file")

    args = parser.parse_args()

    command = args.command

    del args.command

    commands[command](**vars(args))

if __name__ == '__main__':
    main()