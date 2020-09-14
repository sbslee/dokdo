import pandas as pd
import seaborn as sns
import common

@common.visualization
def plot_alpha_rarefaction(**kwargs):

    if not kwargs['method']:
        kwargs['method'] = 'observed'

    methods = {'observed': ['observed_features', 'Observed ASVs'],
               'shannon': ['shannon', 'Shannon Index'],
               'faith': ['faith_pd', "Faith's phylogenetic diversity"]}

    data_file = kwargs['zip_dir'] + '/data/{}.csv'.format(methods[kwargs['method']][0])

    df = pd.read_csv(data_file, index_col=0, sep=',')

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

    sns.lineplot(x='depth', y='ASV', data=data, hue=kwargs['color'],
                 err_style='bars', sort=False, dashes=True, ax=kwargs['ax'])

    kwargs['ax'].set_xlabel('Sequencing depth')
    kwargs['ax'].set_ylabel(methods[kwargs['method']][1])
