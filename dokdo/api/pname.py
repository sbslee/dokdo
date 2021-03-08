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
