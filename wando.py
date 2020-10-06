import os
import argparse

import nbformat as nbf

from create_report import create_report

from compute_table_stat import compute_table_stat

def create_report_(**kwargs):
    create_report(**kwargs)

def pipeline_fastq2asv(**kwargs):
    with open('qsubme-pipeline-fastq2asv.sh', 'w') as f:
        f.write("#$ -S /bin/sh" + '\n')
        f.write("#$ -cwd" + '\n')
        f.write('\n')
        f.write("X_PDP={}".format(kwargs['x_pdp']) + '\n')
        f.write("P_TRUNC_LEN_F={}".format(kwargs['p_trunc_len_f']) + '\n')
        f.write("P_TRUNC_LEN_R={}".format(kwargs['p_trunc_len_r']) + '\n')
        f.write("P_TRIM_LEFT_F={}".format(kwargs['p_trim_left_f']) + '\n')
        f.write("P_TRIM_LEFT_R={}".format(kwargs['p_trim_left_r']) + '\n')
        f.write("P_N_THREADS={}".format(kwargs['p_n_threads']) + '\n')
        f.write('\n')
        f.write(f"sh $X_PDP/pipeline_fastq2asv.sh \\" + '\n')
        f.write("$P_TRUNC_LEN_F \\" + '\n')
        f.write("$P_TRUNC_LEN_R \\" + '\n')
        f.write("$P_TRIM_LEFT_F \\" + '\n')
        f.write("$P_TRIM_LEFT_R \\" + '\n')
        f.write("$P_N_THREADS" + '\n')

def pipeline_analyze(**kwargs):
    with open('qsubme-analyze.sh', 'w') as f:
        f.write("#!/bin/bash" + '\n')
        f.write("#$ -cwd" + '\n')
        f.write('\n')
        f.write(f"sh {os.path.dirname(os.path.realpath(__file__))}/pipeline_analyze.sh {kwargs['p_classifier']}")

def main():
    commands = {'create-report': create_report_,
                'pipeline-fastq2asv': pipeline_fastq2asv,
                'pipeline-analyze': pipeline_analyze,
                'compute-table-stat': compute_table_stat}

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