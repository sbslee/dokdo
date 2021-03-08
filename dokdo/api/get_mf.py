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
    """
    if isinstance(metadata, str):
        mf = Metadata.load(metadata).to_dataframe()
    elif isinstance(metadata, Metadata):
        mf = metadata.to_dataframe()
    else:
        raise TypeError(f"Incorrect metadata type: {type(metadata)}")
    return mf
