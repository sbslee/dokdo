def pname(name):
    """Return a prettified taxon name.

    Parameters
    ----------
    name : str
        Taxon name.
 
    Returns
    -------
    str
        Prettified name.

    Examples
    --------
    Below are some examples.

    >>> a = 'd__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces;s__Schaalia_radingae'
    >>> b = 'Unassigned;__;__;__;__;__;__'
    >>> c = 'd__Bacteria;__;__;__;__;__;__'
    >>> d = 'd__Bacteria;p__Acidobacteriota;c__Acidobacteriae;o__Bryobacterales;f__Bryobacteraceae;g__Bryobacter;__'
    >>>
    >>> print(dokdo.pname(a))
    s__Schaalia_radingae
    >>> print(dokdo.pname(b))
    Unassigned
    >>> print(dokdo.pname(c))
    d__Bacteria
    >>> print(dokdo.pname(d))
    g__Bryobacter
    """
    ranks = list(reversed(name.split(';')))
    for i, rank in enumerate(ranks):
        if rank in ['Others', 'Unassigned']:
            return rank
        if rank == '__':
            continue
        if rank.split('__')[1] is '':
            return ranks[i+1] + ';' + rank
        return rank
