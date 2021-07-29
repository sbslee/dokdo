# Import standard libraries.
import math
import tempfile
import warnings
import numbers

# Import external libraries.
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
import skbio as sb
from skbio.stats.ordination import OrdinationResults
from scipy import stats
from scipy.spatial.distance import euclidean
from skbio.stats.composition import clr

# Import QIIME 2 libraries.
import qiime2
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2 import Visualization
import dokdo

def _pretty_taxa(s):
    """Returns pretty taxa name."""
    if isinstance(s, matplotlib.text.Text):
        s = s.get_text()
    ranks = list(reversed(s.split(';')))

    for i, rank in enumerate(ranks):
        if rank in ['Others', 'Unassigned']:
            return rank

        if rank == '__':
            continue

        if rank.split('__')[1] is '':
            continue

        if 'uncultured' in rank:
            continue

        # The species name can be sometimes tricky to parse because it could
        # be full (e.g. Helicobacter pylori) or partial (e.g. pylori). In the
        # latter case, I will borrow the genus name (e.g. Helicobacter) to
        # form the full species name.
        if 's__' in rank:
            rank = rank.split('__')[1]

            if len(rank.split('_')) == 1:
                genus = ranks[i+1].split('__')[1].split('_')[0]
                species = rank.split('_')[0]
                rank = f'{genus} {species}'
            else:
                rank = rank.replace('_', ' ')

        if '__' in rank:
            rank = rank.split('__')[1]

        return rank

def _artist(
    ax, title=None, title_fontsize=None, xlabel=None, xlabel_fontsize=None,
    ylabel=None, ylabel_fontsize=None, zlabel=None, zlabel_fontsize=None,
    xticks=None, yticks=None, xticklabels=None, xticklabels_fontsize=None,
    yticklabels=None, yticklabels_fontsize=None, xrot=None, xha=None,
    xmin=None, xmax=None, ymin=None, ymax=None, xlog=False, ylog=False,
    hide_xtexts=False, hide_ytexts=False, hide_ztexts=False,
    hide_xlabel=False, hide_ylabel=False, hide_zlabel=False,
    hide_xticks=False, hide_yticks=False, hide_zticks=False,
    hide_xticklabels=False, hide_yticklabels=False, hide_zticklabels=False,
    show_legend=False, legend_loc='best', legend_ncol=1,
    legend_labels=None, remove_duplicates=False, legend_only=False,
    legend_fontsize=None, legend_markerscale=None, legend_lw=None,
    legend_title=None, plot_method=None
):
    """
    This method controls various properties of a figure.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object to draw the plot onto.
    title : str, optional
        Sets the figure title.
    title_fontsize : float or str, optional
        Sets the title font size.
    xlabel : str, optional
        Set the x-axis label.
    xlabel_fontsize : float or str, optional
        Sets the x-axis label font size.
    ylabel : str, optional
        Set the y-axis label.
    ylabel_fontsize : float or str, optional
        Sets the y-axis label font size.
    zlabel : str, optional
        Set the z-axis label.
    zlabel_fontsize : float or str, optional
        Sets the z-axis label font size.
    xticks : list, optional
        Positions of x-axis ticks.
    yticks : list, optional
        Positions of y-axis ticks.
    xticklabels : list, optional
        Tick labels for the x-axis.
    xticklabels_fontsize : float or str, optional
        Font size for the x-axis tick labels.
    yticklabels : list, optional
        Tick labels for the y-axis.
    yticklabels_fontsize : float or str, optional
        Font size for the y-axis tick labels.
    xrot : float, optional
        Rotation degree of tick labels for the x-axis.
    xha : str, optional
        Horizontal alignment of tick labels for the x-axis.
    xmin : float, optional
        Minimum value for the x-axis.
    xmax : float, optional
        Maximum value for the x-axis.
    ymin : float, optional
        Minimum value for the y-axis.
    ymax : float, optional
        Maximum value for the x-axis.
    xlog : bool, default: False
        Draw the x-axis in log scale.
    ylog : bool, default: False
        Draw the y-axis in log scale.
    hide_xtexts : bool, default: False
        Hides all the x-axis texts.
    hide_ytexts : bool, default: False
        Hides all the y-axis texts.
    hide_ztexts : bool, default: False
        Hides all the z-axis texts.
    hide_xlabel : bool, default: False
        Hides the x-axis label.
    hide_ylabel : bool, default: False
        Hides the y-axis label.
    hide_zlabel : bool, default: False
        Hides the z-axis label.
    hide_xticks : bool, default: False
        Hides ticks and tick labels for the x-axis.
    hide_yticks : bool, default: False
        Hides ticks and tick labels for the y-axis.
    hide_zticks : bool, default: False
        Hides ticks and tick labels for the z-axis.
    hide_xticklabels : bool, default: False
        Hides tick labels for the x-axis.
    hide_yticklabels : bool, default: False
        Hides tick labels for the y-axis.
    hide_zticklabels : bool, default: False
        Hides tick labels for the z-axis.
    show_legend : bool, default: False
        Show the figure legend.
    legend_loc : str, default: 'best'
        Legend location specified as in matplotlib.pyplot.legend.
    legend_ncol : int, default: 1
        Number of columns that the legend has.
    legend_only : bool, default: False
        Clear the figure and display the legend only.
    legend_fontsize : float or str, optional
        Sets the legend font size.
    legend_markerscale : float, optional
        Relative size of legend markers compared with the original.
    legend_lw : float, optional
        Width of the lines in the legend.
    legend_title: str, optional
        Legend title.
    plot_method : str, optional
        Name of the plotting method. This argument is internally used for
        the `alpha_rarefaction_plot` method. Not to be used by users.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    Notes
    -----
    Font size can be specified by provding a number or a string as defined in:
    {'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'}.
    """
    if isinstance(title, str):
        ax.set_title(title, fontsize=title_fontsize)

    if isinstance(xlabel, str):
        ax.set_xlabel(xlabel, fontsize=xlabel_fontsize)

    if isinstance(ylabel, str):
        ax.set_ylabel(ylabel, fontsize=ylabel_fontsize)

    if isinstance(zlabel, str):
        ax.set_zlabel(zlabel, fontsize=zlabel_fontsize)

    if isinstance(xticks, list):
        ax.set_xticks(xticks)

    if isinstance(yticks, list):
        ax.set_yticks(yticks)

    if isinstance(xticklabels, list):
        a = len(ax.get_xticklabels())
        b = len(xticklabels)
        if a != b:
            raise ValueError(f"Expected {a} items, but found {b}")
        ax.set_xticklabels(xticklabels)

    if xticklabels_fontsize is not None:
        ax.tick_params(axis='x', which='major', labelsize=xticklabels_fontsize)

    if isinstance(yticklabels, list):
        a = len(ax.get_yticklabels())
        b = len(yticklabels)
        if a != b:
            raise ValueError(f"Expected {a} items, but found {b}")
        ax.set_yticklabels(yticklabels)

    if yticklabels_fontsize is not None:
        ax.tick_params(axis='y', which='major', labelsize=yticklabels_fontsize)

    if isinstance(xrot, numbers.Number):
        ax.set_xticklabels(ax.get_xticklabels(), rotation=xrot)

    if isinstance(xha, str):
        ax.set_xticklabels(ax.get_xticklabels(), ha=xha)

    ax.set_xlim(left=xmin, right=xmax)
    ax.set_ylim(bottom=ymin, top=ymax)

    if xlog:
        ax.set_xscale('log')

    if ylog:
        ax.set_yscale('log')

    if hide_xtexts:
        ax.set_xlabel('')
        ax.set_xticklabels([])

    if hide_ytexts:
        ax.set_ylabel('')
        ax.set_yticklabels([])

    if hide_ztexts:
        ax.set_zlabel('')
        ax.set_zticklabels([])

    if hide_xlabel:
        ax.set_xlabel('')

    if hide_ylabel:
        ax.set_ylabel('')

    if hide_zlabel:
        ax.set_zlabel('')

    if hide_xticks:
        ax.set_xticks([])

    if hide_yticks:
        ax.set_yticks([])

    if hide_zticks:
        ax.set_zticks([])

    if hide_xticklabels:
        ax.set_xticklabels([])

    if hide_yticklabels:
        ax.set_yticklabels([])

    # Control the figure legend.
    h, l = ax.get_legend_handles_labels()

    if legend_labels:
        a = len(legend_labels)
        b = len(l)
        if a != b:
            m = f"Expected {b} legend labels, received {a}"
            raise ValueError(m)
        l = legend_labels

    if remove_duplicates:
        if h:
            n = int(len(h) / 2)
            h, l = h[:n], l[:n]

    def _display_legend():
        leg = ax.legend(h, l, loc=legend_loc, ncol=legend_ncol,
            fontsize=legend_fontsize, markerscale=legend_markerscale,
            title=legend_title, title_fontsize=legend_fontsize)

        if plot_method == 'alpha_rarefaction_plot':
            i = 1
        else:
            i = 0

        if legend_lw is not None:
            for lh in leg.legendHandles[i:]:
                lh.set_linewidth(legend_lw)

    if legend_only:
        # The order matters.
        ax.clear()
        _display_legend()
        ax.axis('off')
    elif show_legend:
        if h:
            _display_legend()
        else:
            warnings.warn("No handles with labels found to put in legend.")
    else:
        if ax.get_legend():
            ax.get_legend().remove()
        else:
            pass

    return ax

def export(input, temp_dir):
    """
    Export QIIME 2 data as files to a temporary directory.

    This method will automatically detect the type of the input file or
    object from QIIME 2 and then export the underlying data as files to the
    specified temporary directory.

    Parameters
    ----------
    input : str, qiime2.Artifact, or qiime2.Visualization
        QIIME 2 file or object.
    temp_dir : str
        Temporary directory.
    """
    if isinstance(input, qiime2.Artifact):
        fn = f'{temp_dir}/temp.qza'
        input.save(fn)
        input = fn
        Artifact.load(input).export_data(temp_dir)
    elif isinstance(input, qiime2.Visualization):
        fn = f'{temp_dir}/temp.qzv'
        input.save(fn)
        input = fn
        Visualization.load(input).export_data(temp_dir)
    elif isinstance(input, str) and input.endswith('.qza'):
        Artifact.load(input).export_data(temp_dir)
    elif isinstance(input, str) and input.endswith('.qzv'):
        Visualization.load(input).export_data(temp_dir)
    else:
        raise TypeError(f'Incorrect input type detected: {type(input)}.')

def get_mf(metadata):
    """
    Convert a Metadata file or object to a dataframe.

    This method automatically detects the type of input metadata and
    then converts it to a :class:`pandas.DataFrame` object.

    Parameters
    ----------
    metadata : str or qiime2.Metadata
        Metadata file or object.

    Returns
    -------
    pandas.DataFrame
        DataFrame object containing metadata.

    Examples
    --------
    This is a simple example.

    >>> mf = dokdo.get_mf('/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv')
    >>> mf.head()
              barcode-sequence  body-site  ...  reported-antibiotic-usage  days-since-experiment-start
    sample-id                              ...
    L1S8          AGCTGACTAGTC        gut  ...                        Yes                          0.0
    L1S57         ACACACTATGGC        gut  ...                         No                         84.0
    L1S76         ACTACGTGTGGT        gut  ...                         No                        112.0
    L1S105        AGTGCGATGCGT        gut  ...                         No                        140.0
    L2S155        ACGATGCGACCA  left palm  ...                         No                         84.0
    """
    if isinstance(metadata, str):
        mf = Metadata.load(metadata).to_dataframe()
    elif isinstance(metadata, Metadata):
        mf = metadata.to_dataframe()
    else:
        raise TypeError(f"Incorrect metadata type: {type(metadata)}")
    return mf

def pname(name, levels=None):
    """
    Return a prettified taxon name.

    Parameters
    ----------
    name : str
        Taxon name.

    Returns
    -------
    str
        Prettified taxon name.

    Examples
    --------

    >>> import dokdo
    >>> dokdo.pname('d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces;s__Schaalia_radingae')
    's__Schaalia_radingae'
    >>> dokdo.pname('Unassigned;__;__;__;__;__;__')
    'Unassigned'
    >>> dokdo.pname('d__Bacteria;__;__;__;__;__;__')
    'd__Bacteria'
    >>> dokdo.pname('d__Bacteria;p__Acidobacteriota;c__Acidobacteriae;o__Bryobacterales;f__Bryobacteraceae;g__Bryobacter;__')
    'g__Bryobacter'
    >>> dokdo.pname('d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces;s__Schaalia_radingae', levels=[6,7])
    'g__Actinomyces;s__Schaalia_radingae'
    """
    if levels is None:
        ranks = list(reversed(name.split(';')))
        for i, rank in enumerate(ranks):
            if rank in ['Others', 'Unassigned']:
                return rank
            if rank == '__':
                continue
            if rank.split('__')[1] is '':
                return ranks[i+1] + ';' + rank
            return rank
    else:
        ranks = name.split(';')
        if 'Others' in ranks:
            return 'Others'
        if 'Unassigned' in ranks:
            return 'Unassigned'
        return ';'.join([ranks[x-1] for x in levels])
