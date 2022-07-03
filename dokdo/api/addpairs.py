import pandas as pd
import seaborn as sns
import numpy as np

def addpairs(
    taxon, csv_file, subject, category, groups,
    ax=None, figsize=None, width=0.8, **kwargs
):
    """
    Add paired lines in a plot created by `dokdo.taxa_abundance_box_plot`.

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
    groups : list
        Groups in the category column.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).
    width : float, default: 0.8
        Width of all the elements for one level of the grouping variable.
    kwargs : other keyword arguments
        All other keyword arguments are passed to :meth:`seaborn.lineplot`.

    Returns
    -------
    matplotlib.axes.Axes
        Axes object with the plot drawn onto it.

    See Also
    --------
    dokdo.api.taxa_abundance_box_plot

    Examples
    --------
    Below is a simple example where we pretend we only have the samples
    shown below and they are from a single subject. We are interested in
    comparing the relative abundance of the phylum Preteobacteria between
    the left palm and right palm. We also want to perofrm the comparison
    in the context of ``days-since-experiment-start`` (i.e. paired
    comparison).

    .. code:: python3

        import dokdo
        from qiime2 import Metadata
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()

        qzv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/taxa-bar-plots.qzv'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'

        sample_names = ['L2S240', 'L3S242', 'L2S155', 'L4S63', 'L2S175', 'L3S313', 'L2S204', 'L4S112', 'L2S222', 'L4S137']

        metadata = Metadata.load(metadata_file)
        metadata = metadata.filter_ids(sample_names)
        mf = dokdo.get_mf(metadata)
        mf = mf[['body-site', 'days-since-experiment-start']]

        ax = dokdo.taxa_abundance_box_plot(
            qzv_file,
            level=2,
            hue='body-site',
            taxa_names=['k__Bacteria;p__Proteobacteria'],
            show_others=False,
            figsize=(8, 7),
            sample_names=sample_names,
            pretty_taxa=True,
            include_samples={'body-site': ['left palm', 'right palm']},
            csv_file='addpairs.csv'
        )

        dokdo.addpairs(
            'k__Bacteria;p__Proteobacteria',
            'addpairs.csv',
            'days-since-experiment-start',
            'body-site',
            ['left palm', 'right palm'],
            ax=ax
        )

        plt.tight_layout()

    .. image:: images/addpairs-1.png

    Note that the method suppors more than two groups.

    .. code:: python3

        metadata = Metadata.load(metadata_file)
        mf = dokdo.get_mf(metadata)
        mf = mf[mf['subject'] == 'subject-1']
        mf = mf[['body-site', 'days-since-experiment-start']]
        mf = mf.drop_duplicates()

        ax = dokdo.taxa_abundance_box_plot(
            qzv_file,
            metadata=Metadata(mf),
            level=2,
            hue='body-site',
            taxa_names=['k__Bacteria;p__Proteobacteria'],
            show_others=False,
            figsize=(8, 7),
            pretty_taxa=True,
            csv_file='addpairs.csv'
        )

        dokdo.addpairs(
            'k__Bacteria;p__Proteobacteria',
            'addpairs.csv',
            'days-since-experiment-start',
            'body-site',
            ['gut', 'left palm', 'right palm', 'tongue'],
            ax=ax
        )

        plt.tight_layout()

    .. image:: images/addpairs-2.png
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    df = pd.read_csv(csv_file)
    df = df[[subject, category, taxon]]
    df = df.rename(columns={taxon: 'perc'})

    n = len(groups)
    half = width / n / 2
    positions = np.linspace(0, width - half * 2, num=n) - half * (n-1)

    d = dict(zip(groups, positions))
    df[category] = df[category].map(d)
    df = df.pivot(index=category, columns=subject, values='perc')
    ax = sns.lineplot(data=df, dashes=False, legend=False, ax=ax, **kwargs)
    ax.set_xlabel('')

    return ax
