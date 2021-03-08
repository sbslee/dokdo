import pandas as pd
from scipy import stats
from .num2sig import num2sig

def wilcoxon(taxon, csv_file, subject, category, group1, group2, ann=False):
    """Compute the p-value from the Wilcoxon Signed-rank test.

    This method tests the null hypothesis that two related paired
    samples come from the same distribution for a given taxon
    using the `scipy.stats.wilcoxon()` method.

    Note that one of the inputs for this method is a .csv file
    from the `dokdo.taxa_abundance_box_plot()` method which
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
    """
    df = pd.read_csv(csv_file)
    df = df.sort_values([subject, category])
    g1 = df[df[category] == group1]
    g2 = df[df[category] == group2]
    p = stats.wilcoxon(g1[taxon], g2[taxon])[1]
    if ann:
        p = num2sig(p)
    return p
