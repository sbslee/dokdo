import pandas as pd
from scipy import stats
from .num2sig import num2sig

def wilcoxon(
    taxon, csv_file, subject, category,
    group1, group2, ann=False
):
    """
    Compute the p-value from the Wilcoxon Signed-rank test.

    This method tests the null hypothesis that two related paired
    samples come from the same distribution for a given taxon
    using the :meth:`scipy.stats.wilcoxon()` method.

    Note that one of the inputs for this method is a .csv file
    from the :meth:`dokdo.taxa_abundance_box_plot()` method which
    contains the relvant data (e.g. relative abundance).

    Parameters
    ----------
    taxon : str
        Target taxon name.
    csv_file : str
        Path to the .csv file.
    subject : str
        Column name to indicate pair information.
    category : str
        Column name to be tested.
    group1 : str
        First group in the category column.
    group2 : str
        Second group in the category column.
    ann : bool, default: False
        If True, return a signifiacne annotation instead of a p-value.
        See `dokdo.num2sig` for how signifiacne levels are defined.

    Returns
    -------
    float or str
        P-value or signifiance annotatiom.

    Examples
    --------
    Below is a simple example where we pretend we only have the samples
    shown below and they are from a single subject. We are interested in
    comparing the relative abundance of the phylum Preteobacteria between
    the left palm and right palm. We also want to perofrm the comparison in
    the context of ``days-since-experiment-start`` (i.e. paired comparison).

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

        taxon = 'k__Bacteria;p__Proteobacteria'
        csv_file = 'wilcoxon.csv'

        ax = dokdo.taxa_abundance_box_plot(
            qzv_file,
            level=2,
            hue='body-site',
            taxa_names=[taxon],
            show_others=False,
            figsize=(6, 6),
            sample_names=sample_names,
            include_samples={'body-site': ['left palm', 'right palm']},
            csv_file=csv_file,
            pretty_taxa=True
        )

        dokdo.addpairs(
            taxon,
            csv_file,
            'days-since-experiment-start',
            'body-site', ['left palm', 'right palm'],
            ax=ax
        )

        plt.tight_layout()

    .. image:: images/wilcoxon.png

    .. code:: python3

        p_value = dokdo.wilcoxon(
            taxon,
            csv_file,
            'days-since-experiment-start',
            'body-site',
            'left palm',
            'right palm'
        )
        print(f'The p-value is {p_value:.6f}')
        # Will print: The p-value is 0.062500
    """
    df = pd.read_csv(csv_file)
    df = df.sort_values([subject, category])
    g1 = df[df[category] == group1]
    g2 = df[df[category] == group2]
    p = stats.wilcoxon(g1[taxon], g2[taxon])[1]
    if ann:
        p = num2sig(p)
    return p
