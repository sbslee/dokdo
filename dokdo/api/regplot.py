import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def regplot(
    taxon, csv_file, subject, category, group1, group2, label=None, ax=None,
    figsize=None
):
    """
    Plot relative abundance data and a linear regression model fit from
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
        Label to use for the legend.
    ax : matplotlib.axes.Axes, optional
        Axes object to draw the plot onto, otherwise uses the current Axes.
    figsize : tuple, optional
        Width, height in inches. Format: (float, float).

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

        sample_names = ['L2S240', 'L3S242', 'L2S155', 'L4S63', 'L2S175', 'L3S313', 'L2S204', 'L4S112', 'L2S222', 'L4S137']

        qzv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/taxa-bar-plots.qzv'
        metadata_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/sample-metadata.tsv'

        metadata = Metadata.load(metadata_file)
        metadata = metadata.filter_ids(sample_names)
        mf = dokdo.get_mf(metadata)
        mf = mf[['body-site', 'days-since-experiment-start']]

        dokdo.taxa_abundance_box_plot(
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

        plt.tight_layout()

    .. image:: images/regplot-1.png

    Finally, run the ``dokdo.regplot`` method.

    .. code:: python3

        dokdo.regplot('k__Bacteria;p__Proteobacteria',
                      'addpairs.csv',
                      'days-since-experiment-start',
                      'body-site',
                      'left palm',
                      'right palm',
                      figsize=(8, 7))

        plt.tight_layout()

    .. image:: images/regplot-2.png
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

    ax.set_xlabel(group1)
    ax.set_ylabel(group2)

    return ax
