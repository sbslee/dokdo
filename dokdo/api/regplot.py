import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from .common import _artist

def regplot(
    taxon, csv_file, subject, category, group1, group2,
    label=None, ax=None, figsize=None, artist_kwargs=None
):
    """Plot relative abundance data and a linear regression model fit from
    paired samples for the given taxon.

    Parameters
    ----------
    taxon : str
        Target taxon name.
    csv_file : str
        Path to the .csv file from the `taxa_abundance_box_plot` method.
    subject : str
        Column name to indicate pair information.
    category : str
        Column name to be studied.
    group1 : str
        First group in the category column.
    group2 : str
        Second group in the category column.
    label : str
        Label to use in a legend.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    taxa_abundance_box_plot

    Examples
    --------
    Below is a simple example where we pretend we only have the samples
    shown below and they are from a single subject. We are interested in
    comparing the relative abundance of the phylum Preteobacteria between
    the left palm and right palm. We also want to perofrm the comparison
    in the context of ``days-since-experiment-start`` (i.e. paired
    comparison).

    .. plot::
        :context: close-figs

        >>> metadata = Metadata.load(f'{data_dir}/moving-pictures-tutorial/sample-metadata.tsv')
        >>> sample_names = ['L2S240', 'L3S242', 'L2S155', 'L4S63', 'L2S175', 'L3S313', 'L2S204', 'L4S112', 'L2S222', 'L4S137']
        >>> metadata = metadata.filter_ids(sample_names)
        >>> mf = dokdo.get_mf(metadata)
        >>> mf = mf[['body-site', 'days-since-experiment-start']]
        >>> mf.sort_values(['days-since-experiment-start', 'body-site'])
                    body-site  days-since-experiment-start
        sample-id
        L2S240      left palm                          0.0
        L3S242     right palm                          0.0
        L2S155      left palm                         84.0
        L4S63      right palm                         84.0
        L2S175      left palm                        112.0
        L3S313     right palm                        112.0
        L2S204      left palm                        140.0
        L4S112     right palm                        140.0
        L2S222      left palm                        168.0
        L4S137     right palm                        168.0

    Next, we will run the ``dokdo.taxa_abundance_box_plot`` method to create the input file for the ``dokdo.regplot`` method.

    .. plot::
        :context: close-figs

        >>> qzv_file = f'{data_dir}/moving-pictures-tutorial/taxa-bar-plots.qzv'
        >>> dokdo.taxa_abundance_box_plot(qzv_file,
        ...                               level=2,
        ...                               hue='body-site',
        ...                               taxa_names=['k__Bacteria;p__Proteobacteria'],
        ...                               show_others=False,
        ...                               figsize=(6, 6),
        ...                               sample_names=sample_names,
        ...                               add_datapoints=True,
        ...                               include_samples={'body-site': ['left palm', 'right palm']},
        ...                               csv_file='addpairs.csv',
        ...                               artist_kwargs=dict(show_legend=True, ymax=70))
        >>> plt.tight_layout()

    Finally, run the ``dokdo.regplot`` method.

    .. plot::
        :context: close-figs

        >>> dokdo.regplot('k__Bacteria;p__Proteobacteria',
        ...               'addpairs.csv',
        ...               'days-since-experiment-start',
        ...               'body-site',
        ...               'left palm',
        ...               'right palm')
        >>> plt.tight_layout()
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    df = pd.read_csv(csv_file)
    df = df.sort_values([subject, category])
    g1 = df[df[category] == group1]
    g2 = df[df[category] == group2]
    df = pd.DataFrame({group1: g1[taxon].to_list(),
                       group2: g2[taxon].to_list()})

    if label is None:
        _label = taxon
    else:
        _label = label

    sns.regplot(data=df, x=group1, y=group2, ax=ax, label=_label)

    if artist_kwargs is None:
        artist_kwargs = {}

    artist_kwargs = {"xlabel": group1,
                     "ylabel": group2,
                     **artist_kwargs}

    ax = _artist(ax, **artist_kwargs)

    return ax
