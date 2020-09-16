import os
import zipfile
import shutil
import common
import pandas as pd
import numpy as np

def compute_table_stat(**kwargs):
    kwargs = {**common.KWARGS, **kwargs}

    if kwargs['i_zip']:
        with zipfile.ZipFile(kwargs['i_zip'], 'r') as zip_file:
            zip_file.extractall()
            kwargs['zip_dir'] = zip_file.namelist()[0].split('/')[0]

    a = kwargs['zip_dir'] + f'/data/sample-frequency-detail.csv'

    df = pd.read_csv(a, sep=',', header=None)

    a = df.iloc[:, 1].astype(int)

    if kwargs['p_stat'] in [None, 'minimum']:
        print(a.min())
    elif kwargs['p_stat'] == 'maximum':
        print(a.max())
    elif kwargs['p_stat'] == 'median':
        print(a.median())
    elif kwargs['p_stat'] == 'mean':
        print(a.mean())
    else:
        raise ValueError("Incorrect --p-stat detected")

    if not kwargs['keep']:
        shutil.rmtree(kwargs['zip_dir'])