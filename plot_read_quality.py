import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import common

@common.visualization
def plot_read_quality(**kwargs):

    if kwargs['p_forward']:
        a = 'forward'
    else:
        a = 'reverse'

    a = kwargs['zip_dir'] + f'/data/{a}-seven-number-summaries.tsv'

    df = pd.read_csv(a, index_col=0, sep='\t', skiprows=[1])
    df = pd.melt(df.reset_index(), id_vars=['index'])
    df['variable'] = df['variable'].astype('int64')

    sns.boxplot(x='variable', y='value', data=df, ax=kwargs['ax'],
                fliersize=0, color='white', whiskerprops={'linestyle': ':'})

    kwargs['ax'].set_ylim([0, 45])
    kwargs['ax'].set_xlabel('Sequence base')
    kwargs['ax'].set_ylabel('Quality score')

    a = np.arange(df['variable'].min(), df['variable'].max(), 20)
    kwargs['ax'].set_xticks(a)
    kwargs['ax'].set_xticklabels(a)