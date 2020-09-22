from tempfile import TemporaryDirectory

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from qiime2 import Metadata
from qiime2 import Visualization

def ancom_volcano_plot(ancom, ax=None):
    """
    This method creates an ANCOM volcano plot.

    Example:
        api.ancom_volcano_plot('ancom-Site.qzv')
    """
    t = TemporaryDirectory()
    Visualization.load(ancom).export_data(t.name)
    df = pd.read_table(f'{t.name}/data.tsv')
    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 10))
    ax.scatter(df.clr, df.W, s=80, c='black', alpha=0.5)
    ax.set_xlabel('clr')
    ax.set_ylabel('W')



def alpha_diversity_plot(significance, where, ax=None):
    """
    This method creates an alpha diversity plot.

    Example:
        api.alpha_diversity_plot('shannon_group-significance.qzv', 'Site')
    """
    t = TemporaryDirectory()
    Visualization.load(significance).export_data(t.name)
    df = Metadata.load(f'{t.name}/metadata.tsv').to_dataframe()
    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 10))
    metric = df.columns[-1]
    boxprops = dict(color='white', edgecolor='black')
    sns.boxplot(x=where, y=metric, data=df, ax=ax, boxprops=boxprops)
