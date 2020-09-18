import os
import subprocess

from qiime2.plugins import feature_table
from qiime2 import Artifact
import pandas as pd

def action_collapse(i_table, i_tax, o_prefix, **kwargs):
    for i in range(1, 8):
        n = str(i)
        fn = f'{o_prefix}-level-{n}.qza'

        command = ['qiime', 'taxa', 'collapse',
                   '--i-table', 'table.qza',
                   '--i-taxonomy', 'taxonomy.qza',
                   '--p-level', n,
                   '--o-collapsed-table', fn]

        subprocess.run(command)

        table = Artifact.load(fn)
        df = table.view(pd.DataFrame)
        df = df.sort_index().T
        df.to_csv(fn.replace('qza', 'csv'))
        os.remove(fn)