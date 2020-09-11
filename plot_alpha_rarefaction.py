import os
import zipfile
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_alpha_rarefaction(qzv_file,
                           color,
                           method='observed',
                           ax=None,
                           show=False,
                           figsize=None,
                           legend_fontsize=8,
                           tick_fontsize=8,
                           label_fontsize=8,
                           keep=False,
                           output=None):

    methods = {'observed': ['observed_features', 'Observed ASVs'],
               'shannon': ['shannon', 'Shannon Index'],
               'faith': ['faith_pd', "Faith's phylogenetic diversity"]}

    with zipfile.ZipFile(qzv_file, 'r') as zip_file:
        zip_file.extractall()
        zip_dir = zip_file.namelist()[0].split('/')[0]
        csv_file = zip_dir + f'/data/{methods[method][0]}.csv'

    df = pd.read_csv(csv_file, index_col=0, sep=',')

    if not keep:
        shutil.rmtree(zip_dir)

    cols = [x for x in df.columns if 'iter' not in x]
    mean = df[cols]
    data = pd.DataFrame(columns=cols)
    depths = [col.split('_')[0] for col in df.columns if 'depth' in col]
    df = df.drop(cols,axis=1)
    df.columns = depths

    for depth in depths:
        mean['ASV'] = df[depth].mean(axis=1)
        mean['depth']= depth.split('-')[-1]
        data = pd.concat([data, mean], sort=True)

    if not ax:
        fig, ax = plt.subplots(figsize=figsize)

    sns.lineplot(x='depth', y='ASV', data=data, hue=color,
                 err_style='bars', sort=False, dashes=True, ax=ax)

    ax.set_xlabel('Sequencing depth', fontsize=label_fontsize)
    ax.set_ylabel(f'{methods[method][1]}', fontsize=label_fontsize)
    ax.tick_params(axis='both', labelsize=tick_fontsize)
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=legend_fontsize)

    if show:
        plt.show()

    if output:
        plt.savefig(output)