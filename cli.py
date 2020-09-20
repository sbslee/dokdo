import argparse
import os

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

def make_manifest(fastq_dir, output):
    files = {}

    for r, d, f in os.walk(fastq_dir):
        for x in f:
            name = x.split('_')[0]

            if '_R1_001.fastq' in x:
                if name not in files:
                    files[name] = ['', '']
                files[name][0] = f'{r}/{x}'
            elif '_R2_001.fastq' in x:
                if name not in files:
                    files[name] = ['', '']
                files[name][1] = f'{r}/{x}'
            else:
                pass

    with open(output, 'w') as f:
        headers = ['sample-id', 'forward-absolute-filepath', 
                   'reverse-absolute-filepath']
        f.write('\t'.join(headers) + '\n')

        for name in sorted(files):
            fields = [name, files[name][0], files[name][1]]
            f.write('\t'.join(fields) + '\n')

def main():
    commands = {
        'collapse': collapse,
        'tax2seq': tax2seq,
        'make_manifest': make_manifest,
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

    make_manifest_parser = subparsers.add_parser(
        'make_manifest',
        help="create manifest file from FASTQ directory",
    )
    make_manifest_parser.add_argument('fastq_dir',
                                 help="path to input FASTQ directory")
    make_manifest_parser.add_argument('output',
                                 help="path to output manifest file")

    args = parser.parse_args()

    command = args.command

    del args.command

    commands[command](**vars(args))

if __name__ == '__main__':
    main()
