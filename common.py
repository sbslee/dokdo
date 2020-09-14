import os
import zipfile
import shutil
import matplotlib.pyplot as plt

TAXA = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

KWARGS = {'figsize': None,
          'keep': False,
          'ax': None,
          'xlabel_fontsize': None,
          'ylabel_fontsize': None,
          'tick_fontsize': None,
          'legend_fontsize': None,
          'legend': True,
          'show': False,
          'output': None,
          'legend_fontsize': None,
          'zip_dir': None,
          'qzv_file': None,
          'color': None,
          'width': 0.8,
          'method': None,
          'sortby': []}

def visualization(func):
    def wrapper(*args, **kwargs):
        kwargs = {**KWARGS, **kwargs}

        with zipfile.ZipFile(kwargs['qzv_file'], 'r') as zip_file:
            zip_file.extractall()
            kwargs['zip_dir'] = zip_file.namelist()[0].split('/')[0]

        if not kwargs['ax']:
            fig, kwargs['ax'] = plt.subplots(figsize=kwargs['figsize'])

        func(*args, **kwargs)

        kwargs['ax'].tick_params(axis='both',
                                labelsize=kwargs['tick_fontsize'])

        kwargs['ax'].xaxis.label.set_fontsize(kwargs['xlabel_fontsize'])

        kwargs['ax'].yaxis.label.set_fontsize(kwargs['ylabel_fontsize'])

        if kwargs['legend']:
            kwargs['ax'].legend(loc='center left',
                                bbox_to_anchor=(1, 0.5),
                                fontsize=kwargs['legend_fontsize'])

        plt.tight_layout()

        if kwargs['output']:
            plt.savefig(kwargs['output'])

        if kwargs['show']:
            plt.show()

        if not kwargs['keep']:
            shutil.rmtree(kwargs['zip_dir'])

    return wrapper