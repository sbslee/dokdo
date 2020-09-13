import os
import zipfile
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import common

def plot_taxa_abundance(qzv_file,
                        sortby=[],
                        show=False,
                        color=None,
                        method='Relative',
                        figsize=(20, 10),
                        legend=True,
                        legend_fontsize=8,
                        xlabel_fontsize=8,
                        ylabel_fontsize=8,
                        tick_fontsize=8,
                        width=0.9,
                        output=None):

    if not color:
        color = 'Kingdom'

    level = common.TAXA.index(color) + 1

    methods = {'Relative': 'Relative frequency',
               'Absolute': 'Frequency'}

    with zipfile.ZipFile(qzv_file, 'r') as zip_file:
        zip_file.extractall()
        zip_dir = zip_file.namelist()[0].split('/')[0]
        csv_file = zip_dir + f'/data/level-{level}.csv'

    df = pd.read_csv(csv_file, index_col=0, sep=',')

    if sortby:
        df = df.sort_values(by=sortby)

    dropped = []

    for column in df.columns:
        if 'Unassigned' in column:
            continue
        elif '__' in column:
            continue
        else:
            dropped.append(column)

    df = df.drop(columns=dropped)

    if method == 'Relative':
        df = df.T
        df = df / df.sum()
        df = df.T

    df = df.loc[:, df.mean().sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=figsize)

    df.plot.bar(stacked=True, legend=False, ax=ax, width=width,
                color=plt.cm.get_cmap('Accent').colors)

    ax.set_xlabel('Samples', fontsize=xlabel_fontsize)
    ax.set_ylabel(f'{methods[method]}', fontsize=ylabel_fontsize)
    ax.tick_params(axis='both', labelsize=tick_fontsize)

    plt.tight_layout()

    if legend:
        ax.legend(loc='center left',
                  bbox_to_anchor=(1, 0.5),
                  fontsize=legend_fontsize)

    if show:
        plt.show()

    if output:
        plt.savefig(output)