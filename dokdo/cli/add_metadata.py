import pandas as pd
from qiime2 import Metadata
import warnings

def add_metadata(metadata_file,
                 columns_file,
                 output_file):
    """Add new metadata columns to an existing sample-metadata file (.tsv).

    Parameters
    ----------

    Path to the output file.
    """
    mf1 = Metadata.load(metadata_file).to_dataframe()
    index_name = mf1.index.name
    dtypes = mf1.dtypes.to_dict()
    mf2 = pd.read_table(columns_file, keep_default_na=False)

    for k, v in dtypes.items():
        if k in mf2.columns:
            if v == 'object':
                mf2[k] = mf2[k].astype(str)
            else:
                mf2[k] = mf2[k].astype(v)

    mf3 = mf1.reset_index().merge(mf2).set_index(index_name)
    mf3 = mf3.reindex(mf1.index)

    a = mf1.shape[0]
    b = mf3.shape[0]

    if a != b:
        message = (f"Final metadata (N={b}) has different number of samples "
                   f"than input metadata (N={a}). Please double check "
                   "whether this was intended.")
        warnings.warn(message)

    if mf3.isnull().values.any():
        warnings.warn("Final metadata contains NaN. Please double check "
                      "whether this was intended.")

    Metadata(mf3).save(output_file)
