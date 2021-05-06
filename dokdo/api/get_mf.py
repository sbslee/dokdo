from qiime2 import Metadata

def get_mf(metadata):
    """Convert a file or object from QIIME 2 metadata to a dataframe.

    This method automatically detects the type of input metadata and
    then converts it to a pandas.DataFrame object.

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
