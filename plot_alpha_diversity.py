import pandas as pd
import seaborn as sns
import common

@common.visualization
def plot_alpha_diversity(**kwargs):

    methods = {'pielou_evenness': "Pielou's evenness index",
               'shannon_entropy': 'Shannon index',
               'faith_pd': "Faith's phylogenetic diversity"}

    df = pd.read_csv(kwargs['zip_dir'] + f'/data/metadata.tsv',
                     index_col=0, sep='\t', skiprows=[1])

    df = df.sort_values(by=kwargs['color'])

    for k in methods:
        if k in df.columns:
            method = k
            break

    sns.boxplot(x=kwargs['color'], y=method, data=df, ax=kwargs['ax'])

    kwargs['ax'].set_xlabel(kwargs['color'])
    kwargs['ax'].set_ylabel(f'{methods[method]}')
