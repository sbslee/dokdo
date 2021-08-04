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
