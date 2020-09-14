import os
import argparse

from plot_alpha_rarefaction import plot_alpha_rarefaction
from plot_taxa_abundance import plot_taxa_abundance
from plot_alpha_diversity import plot_alpha_diversity

def make_manifest(i_path=None, o_path=None, **kwargs):
    files = {}

    for r, d, f in os.walk(i_path):
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

    with open(o_path, 'w') as f:
        headers = ['sample-id', 'forward-absolute-filepath', 
                   'reverse-absolute-filepath']
        f.write('\t'.join(headers) + '\n')

        for name in sorted(files):
            fields = [name, files[name][0], files[name][1]]
            f.write('\t'.join(fields) + '\n')

def merge_metadata(i_paths=None, o_path=None, **kwargs):
    metadata = []

    for i in range(len(i_paths)):
        with open(i_paths[i]) as f:
            for j, line in enumerate(f):
                fields = line.strip().split('\t')
                if i == 0 and j < 2:
                    metadata.append(fields)

                if j < 2:
                    continue

                metadata.append(fields)

    with open(o_path, 'w') as f:
        for fields in metadata:
            f.write('\t'.join(fields) + '\n')

def plot_alpha_rarefaction_(i_path=None,
                            o_path=None,
                            p_color=None, 
                            p_figsize=None,
                            **kwargs):

    plot_alpha_rarefaction(i_path,
                           output=o_path,
                           color=p_color,
                           figsize=p_figsize)

def plot_taxa_abundance_(i_path=None,
                         o_path=None,
                         p_color=None,
                         p_figsize=None,
                         **kwargs):

    plot_taxa_abundance(i_path,
                        output=o_path,
                        color=p_color,
                        figsize=p_figsize)

def plot_alpha_diversity_(**kwargs):

    plot_alpha_diversity(qzv_file=kwargs['i_path'],
                         output=kwargs['o_path'],
                         figsize=kwargs['p_figsize'],
                         color=kwargs['p_color'])

def fastq2asv(i_path, p_trunc_len_f=None, p_trunc_len_r=None, p_trim_left_f=None, p_trim_left_r=None, **kwargs):
    with open('qsubme.sh', 'w') as f:
        f.write("#!/bin/bash" + '\n')
        f.write("#$ -cwd" + '\n')
        f.write('\n')
        f.write(f"sh {os.path.dirname(os.path.realpath(__file__))}/fastq2asv.sh {i_path} {p_trunc_len_f} {p_trunc_len_r} {p_trim_left_f} {p_trim_left_r}")

def main():
    commands = {'make-manifest': make_manifest,
                'merge-metadata': merge_metadata,
                'plot-alpha-rarefaction': plot_alpha_rarefaction_,
                'plot-taxa-abundance': plot_taxa_abundance_,
                'plot-alpha-diversity': plot_alpha_diversity_,
                'fastq2asv': fastq2asv}

    parser = argparse.ArgumentParser()

    parser.add_argument('--command')
    parser.add_argument('--i-path')
    parser.add_argument('--i-paths', action='append')
    parser.add_argument('--p-color')
    parser.add_argument('--p-trunc-len-f')
    parser.add_argument('--p-trunc-len-r')
    parser.add_argument('--p-trim-left-f')
    parser.add_argument('--p-trim-left-r')
    parser.add_argument('--p-figsize', nargs=2, type=int)
    parser.add_argument('--o-path')

    args = parser.parse_args()

    commands[args.command](**vars(args))

if __name__ == '__main__':
    main()