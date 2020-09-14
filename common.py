import os
import zipfile
import shutil
import matplotlib.pyplot as plt

TAXA = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

KWARGS = {'p_figsize': None,
          'keep': False,
          'ax': None,
          'xlabel_fontsize': None,
          'ylabel_fontsize': None,
          'tick_fontsize': None,
          'p_legend_title': None,
          'legend_fontsize': None,
          'legend': True,
          'show': False,
          'output': None,
          'zip_dir': None,
          'i_path': None,
          'p_color': None,
          'width': 0.8,
          'method': None,
          'sortby': [],
          'm_path': None,
          'p_dim': 2,
          'p_elev': None,
          'p_azim': None,
          'p_s': 20}

def visualization(func):
    def wrapper(*args, **kwargs):
        kwargs = {**KWARGS, **kwargs}

        with zipfile.ZipFile(kwargs['i_path'], 'r') as zip_file:
            zip_file.extractall()
            kwargs['zip_dir'] = zip_file.namelist()[0].split('/')[0]

        if not kwargs['ax']:
            if kwargs['p_dim'] == 2:
                fig, kwargs['ax'] = plt.subplots(figsize=kwargs['p_figsize'])
            elif kwargs['p_dim'] == 3:
                fig = plt.figure(figsize=kwargs['p_figsize'])
                kwargs['ax'] = fig.add_subplot(1, 1, 1, projection='3d')
            else:
                raise ValueError("Incorrect dimension detected")

        if kwargs['p_dim'] == 3:
            kwargs['ax'].view_init(elev=kwargs['p_elev'],
                                   azim=kwargs['p_azim'])

        func(*args, **kwargs)

        kwargs['ax'].tick_params(axis='both',
                                labelsize=kwargs['tick_fontsize'])

        kwargs['ax'].xaxis.label.set_fontsize(kwargs['xlabel_fontsize'])

        kwargs['ax'].yaxis.label.set_fontsize(kwargs['ylabel_fontsize'])

        if kwargs['legend']:
            if kwargs['p_legend_title']:
                kwargs['ax'].legend(title=kwargs['p_legend_title'],
                                    loc='center left',
                                    bbox_to_anchor=(1, 0.5),
                                    fontsize=kwargs['legend_fontsize'])

            else:
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