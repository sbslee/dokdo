import os
import argparse

from compute_table_stat import compute_table_stat

def main():
    commands = {'compute-table-stat': compute_table_stat}

    parser = argparse.ArgumentParser()

    parser.add_argument('--command')
    parser.add_argument('--i-path')
    parser.add_argument('--i-paths', action='append')
    parser.add_argument('--i-zip')
    parser.add_argument('--i-fastq')
    parser.add_argument('--i-rep-seqs')
    parser.add_argument('--i-tax')
    parser.add_argument('--i-table')
    parser.add_argument('--p-color')
    parser.add_argument('--p-trunc-len-f', default=0, type=int)
    parser.add_argument('--p-trunc-len-r', default=0, type=int)
    parser.add_argument('--p-trim-left-f', default=0, type=int)
    parser.add_argument('--p-trim-left-r', default=0, type=int)
    parser.add_argument('--p-n-threads', default=1, type=int)
    parser.add_argument('--p-figsize', nargs=2, type=int)
    parser.add_argument('--p-stat')
    parser.add_argument('--p-classifier')
    parser.add_argument('--p-forward', action='store_true')
    parser.add_argument('--p-tax')
    parser.add_argument('--o-path')
    parser.add_argument('--o-prefix')

    args = parser.parse_args()

    args.x_pdp = os.path.dirname(os.path.realpath(__file__))

    commands[args.command](**vars(args))

if __name__ == '__main__':
    main()