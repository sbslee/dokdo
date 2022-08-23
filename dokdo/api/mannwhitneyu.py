import pandas as pd
from scipy import stats
from .num2sig import num2sig

def mannwhitneyu(
    taxon, csv_file, category,
    group1, group2, ann=False
):
    """
    Compute the p-value from the Mann–Whitney U test.

    This method tests the null hypothesis that two independent
    samples come from the same distribution for a given taxon
    using the :meth:`scipy.stats.mannwhitneyu()` method.

    Note that one of the inputs for this method is a .csv file
    from the :meth:`dokdo.taxa_abundance_box_plot()` method which
    contains the relvant data (e.g. relative abundance).

    Parameters
    ----------
    taxon : str
        Target taxon name.
    csv_file : str
        Path to the .csv file.
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
        P-value or signifiance annotation.

    Examples
    --------
    Below is a sample example where we compare the relative abundance of
    the phylum Proteobacteria between the left palm and right palm.
    Before we can calculate the p-value for this comparison, we first
    need to create a .csv file containing the relevant data using the
    ``dokdo.taxa_abundance_box_plot`` method.

    .. code:: python3

        import dokdo
        from qiime2 import Metadata
        import matplotlib.pyplot as plt
        %matplotlib inline
        import seaborn as sns
        sns.set()

        from qiime2 import Metadata
        taxon = 'k__Bacteria;p__Proteobacteria'
        csv_file = 'mannwhitneyu.csv'
        qzv_file = '/Users/sbslee/Desktop/dokdo/data/moving-pictures-tutorial/taxa-bar-plots.qzv'

        ax = dokdo.taxa_abundance_box_plot(
            qzv_file,
            level=2,
            hue='body-site',
            taxa_names=[taxon],
            show_others=False,
            figsize=(8, 7),
            pretty_taxa=True,
            include_samples={'body-site': ['left palm', 'right palm']},
            csv_file=csv_file
        )

        plt.tight_layout()

    .. image:: images/mannwhitneyu.png

    .. code:: python3

        p_value = dokdo.mannwhitneyu(
            taxon,
            csv_file,
            'body-site',
            'left palm',
            'right palm'
        )
        print(f'The p-value is {p_value:.6f}')
        # Will print: The p-value is 0.235243
    """
    df = pd.read_csv(csv_file)
    df = df[[taxon, category]]
    df = df.rename(columns={taxon: 'perc'})
    df = df.sort_values([category])
    g1 = df[df[category] == group1]
    g2 = df[df[category] == group2]
    p = stats.mannwhitneyu(g1['perc'], g2['perc'])[1]
    if ann:
        p = num2sig(p)
    return p
