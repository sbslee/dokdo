import pandas as pd
from scipy import stats
from .num2sig import num2sig

def mannwhitneyu(taxon, csv_file, category, group1, group2, ann=False):
    """Compute the p-value from the Mann–Whitney U test.

    This method tests the null hypothesis that two independent samples
    come from the same distribution for a given taxon.

    Parameters
    ----------
    taxon : str
        Target taxon name.
    csv_file : str
        Path to the .csv file from `dokdo.taxa_abundance_box_plot`.
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
        The p-value or signifiance annotation for the test
        from `scipy.stats.mannwhitneyu`.
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
