from tempfile import TemporaryDirectory

import pandas as pd
import matplotlib.pyplot as plt

from qiime2 import Visualization

def ancom_volcano_plot(ancom, ax=None):
    """
    This method creates an ANCOM volcano plot.

    Example:
        import matplotlib.pyplot as plt
        api.ancom_volcano_plot('ancom-Site.qzv')
        plt.savefig('ancom-Site.pdf')
    """
    t = TemporaryDirectory()
    Visualization.load(ancom).export_data(t.name)
    df = pd.read_table(f'{t.name}/data.tsv')
    if ax is None:
        fig, ax = plt.subplots(figsize=(15, 10))
    ax.scatter(df.clr, df.W, s=80, c='black', alpha=0.5)
    ax.set_xlabel('clr')
    ax.set_ylabel('W')
