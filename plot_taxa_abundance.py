import pandas as pd
import matplotlib.pyplot as plt
import common

@common.visualization
def plot_taxa_abundance(**kwargs):

    methods = {'Relative': 'Relative frequency',
               'Absolute': 'Frequency'}

    if not kwargs['color']:
        kwargs['color'] = 'Kingdom'

    if not kwargs['method']:
        kwargs['method'] = 'Relative'

    level = common.TAXA.index(kwargs['color']) + 1

    df = pd.read_csv(kwargs['zip_dir'] + f'/data/level-{level}.csv',
                     index_col=0, sep=',')

    if kwargs['sortby']:
        df = df.sort_values(by=kwargs['sortby'])

    dropped = []

    for column in df.columns:
        if 'Unassigned' in column:
            continue
        elif '__' in column:
            continue
        else:
            dropped.append(column)

    df = df.drop(columns=dropped)

    if kwargs['method'] == 'Relative':
        df = df.T
        df = df / df.sum()
        df = df.T

    df = df.loc[:, df.mean().sort_values(ascending=False).index]

    df.plot.bar(stacked=True,
                legend=False,
                ax=kwargs['ax'],
                width=kwargs['width'],
                color=plt.cm.get_cmap('Accent').colors)

    kwargs['ax'].set_xlabel('Samples')
    kwargs['ax'].set_ylabel(methods[kwargs['method']])
