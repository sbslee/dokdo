import os
import zipfile
import shutil
import matplotlib.pyplot as plt

TAXA = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

KWARGS = {'i_path': None,
          'i_zip': None,
          'p_figsize': None,
          'p_title': None,
          'p_xlabel_fontsize': None,
          'p_ylabel_fontsize': None,
          'p_tick_fontsize': None,
          'p_legend_title': None,
          'p_legend_fontsize': None,
          'p_legend': True,
          'p_dim': 2,
          'p_elev': None,
          'p_azim': None,
          'p_s': 20,
          'p_color': None,
          'p_forward': True,
          'show': False,
          'output': None,
          'zip_dir': None,
          'width': 0.8,
          'method': None,
          'sortby': [],
          'm_path': None,
          'keep': False,
          'ax': None}

def visualization(func):
    def wrapper(*args, **kwargs):
        kwargs = {**KWARGS, **kwargs}

        if kwargs['i_zip']:
            with zipfile.ZipFile(kwargs['i_zip'], 'r') as zip_file:
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
                                labelsize=kwargs['p_tick_fontsize'])

        kwargs['ax'].xaxis.label.set_fontsize(kwargs['p_xlabel_fontsize'])

        kwargs['ax'].yaxis.label.set_fontsize(kwargs['p_ylabel_fontsize'])

        # Manage the main title.
        kwargs['ax'].title.set_text(kwargs['p_title'])

        # Manage the legend.
        if kwargs['p_legend']:
            handles, labels = kwargs['ax'].get_legend_handles_labels()

            if kwargs['p_legend_title']:
                labels[0] = kwargs['p_legend_title']

                kwargs['ax'].legend(handles, labels,
                                    loc='center left',
                                    bbox_to_anchor=(1, 0.5),
                                    fontsize=kwargs['p_legend_fontsize'])

        # Tighten up the plot.
        plt.tight_layout()

        if kwargs['output']:
            plt.savefig(kwargs['output'])

        if kwargs['show']:
            plt.show()

        if not kwargs['keep']:
            shutil.rmtree(kwargs['zip_dir'])

    return wrapper