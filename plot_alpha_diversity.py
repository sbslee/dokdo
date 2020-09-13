import os
import zipfile
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import common

def plot_alpha_diversity(qzv_file,
                         color,
                         **kwargs):

    kwargs = {**common.KWARGS, **kwargs}

    methods = {'pielou_evenness': "Pielou's evenness index",
               'shannon_entropy': 'Shannon index',
               'faith_pd': "Faith's phylogenetic diversity"}

    with zipfile.ZipFile(qzv_file, 'r') as zip_file:
        zip_file.extractall()
        zip_dir = zip_file.namelist()[0].split('/')[0]
        csv_file = zip_dir + f'/data/metadata.tsv'

    df = pd.read_csv(csv_file, index_col=0, sep='\t', skiprows=[1])
    df = df.sort_values(by=color)

    if not kwargs['keep']:
        shutil.rmtree(zip_dir)

    if not kwargs['ax']:
        fig, ax = plt.subplots(figsize=kwargs['figsize'])

    for k in methods:
        if k in df.columns:
            method = k
            break

    sns.boxplot(x=color, y=method, data=df, ax=ax)

    ax.set_xlabel(color, fontsize=kwargs['xlabel_fontsize'])
    ax.set_ylabel(f'{methods[method]}', fontsize=kwargs['ylabel_fontsize'])
    ax.tick_params(axis='both', labelsize=kwargs['tick_fontsize'])

    plt.tight_layout()

    if kwargs['legend']:
        ax.legend(loc='center left',
                  bbox_to_anchor=(1, 0.5),
                  fontsize=kwargs['legend_fontsize'])

    if kwargs['show']:
        plt.show()

    if kwargs['output']:
        plt.savefig(output)
