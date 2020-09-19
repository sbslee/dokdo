import argparse

import pandas as pd

from qiime2 import Artifact
from q2_types.feature_data import DNAFASTAFormat

from qiime2.plugins import taxa

def collapse(table, taxonomy):
    for i in range(1, 8):
        _ = taxa.methods.collapse(table=Artifact.load(table),
                                  taxonomy=Artifact.load(taxonomy),
                                  level=i)
        _.collapsed_table.view(pd.DataFrame).T.to_csv(f'level-{i}.csv')

def tax2seq(taxonomy, rep_seqs, output):
    tax_df = Artifact.load(taxonomy).view(pd.DataFrame)

    _ = Artifact.load(rep_seqs)
    f = _.view(DNAFASTAFormat).open()

    features = []
    seqs = []

    for i, line in enumerate(f):
        if i % 2 == 0:
            features.append(line.strip().replace('>', ''))
        else:
            seqs.append(line.strip())

    seq_df = pd.DataFrame({'Feature ID': features, 'Sequence': seqs})
    seq_df = seq_df.set_index('Feature ID')

    _ = pd.concat([tax_df, seq_df], axis=1, sort=False)
    _.to_csv(output)

def main():
    commands = {
        'collapse': collapse,
        'tax2seq': tax2seq,
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

    tax2seq_parser = subparsers.add_parser(
        'tax2seq',
        help="return mapping between ASVs and taxonomic classifications",
    )
    tax2seq_parser.add_argument('taxonomy',
                                help="path to input taxonomy.qza file")
    tax2seq_parser.add_argument('rep_seqs',
                                help="path to input rep-seqs.qza file")
    tax2seq_parser.add_argument('output',
                                help="path to output csv file")

    args = parser.parse_args()

    command = args.command

    del args.command

    commands[command](**vars(args))

if __name__ == '__main__':
    main()
