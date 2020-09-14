import pandas as pd
import seaborn as sns
import common

@common.visualization
def plot_read_quality(**kwargs):

    data_file = kwargs['zip_dir'] + '/data/forward-seven-number-summaries.tsv'
    df = pd.read_csv(data_file, index_col=0, sep='\t', skiprows=[1])
    df = pd.melt(df.reset_index(), id_vars=['index'])

    df['variable'] = df['variable'].astype('int64')
    print(df)

#     sns.lineplot(x='variable', y='value', data=df,
#                  err_style='bars', sort=False, dashes=True, ax=kwargs['ax'])

    sns.boxplot(x='variable', y='value', data=df, fliersize=0, color="seagreen", whiskerprops = dict(linestyle='-',linewidth=1.0, color='black'))

    kwargs['ax'].set_ylim([0, 45])

    kwargs['ax'].set_xlabel('Sequence base')
    kwargs['ax'].set_ylabel('Quality score')

    every_nth = 20
    for n, label in enumerate(kwargs['ax'].xaxis.get_ticklabels()):
        if n % every_nth != 0:
            label.set_visible(False)