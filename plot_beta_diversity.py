import pandas as pd
import matplotlib.pyplot as plt
import common

@common.visualization
def plot_beta_diversity(**kwargs):

    data_file = kwargs['zip_dir'] + '/data/ordination.txt'

    df = pd.read_table(data_file,
                       skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8],
                       skipfooter=4,
                       engine='python',
                       header=None,
                       index_col=0)

    df = df.sort_index()

    m_df = pd.read_table(kwargs['m_path'], skiprows=[1], index_col=0)
    
    m_df = m_df.sort_index()

    for c in sorted(m_df[kwargs['p_color']].unique()):
        i = m_df[kwargs['p_color']] == c
        df_list = [df[i].iloc[:, j] for j in range(kwargs['p_dim'])]
        kwargs['ax'].scatter(*df_list, label=c, s=kwargs['p_s'])
