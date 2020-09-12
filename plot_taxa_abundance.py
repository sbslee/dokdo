import os
import zipfile
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import common

def plot_taxa_abundance(qzv_file,
                        color='Kingdom',
                        method='relative',
                        figsize=(20, 10),
                        legend=True,
                        legend_fontsize=8,
                        xlabel_fontsize=8,
                        ylabel_fontsize=10,
                        tick_fontsize=8):

    level = common.TAXA.index(color) + 1

    methods = {'relative': 'Relative frequency',
               'absolute': 'Frequency'}

    with zipfile.ZipFile(qzv_file, 'r') as zip_file:
        zip_file.extractall()
        zip_dir = zip_file.namelist()[0].split('/')[0]
        csv_file = zip_dir + f'/data/level-{level}.csv'

    df = pd.read_csv(csv_file, index_col=0, sep=',')
    df = df.drop(columns=['Site', 'Set'])

    if method == 'relative':
        df = df.T
        df = df / df.sum()
        df = df.T

    fig, ax = plt.subplots(figsize=figsize)
    #df.sort_index().plot(kind='bar', cmap='Accent', stacked=True, legend=False, ax=ax)
    df.sort_index().plot.bar(stacked=True, legend=False, ax=ax)
    ax.set_xlabel('Samples', fontsize=xlabel_fontsize)
    ax.set_ylabel(f'{methods[method]}', fontsize=ylabel_fontsize)

    ax.tick_params(axis='both', labelsize=tick_fontsize)

    plt.tight_layout()

    if legend:
        ax.legend(loc='center left',
                  bbox_to_anchor=(1, 0.5),
                  fontsize=legend_fontsize)