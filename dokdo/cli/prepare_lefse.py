import dokdo
import pandas as pd
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins import feature_table
from qiime2.plugins import taxa

def prepare_lefse(table_file,
                  taxonomy_file,
                  metadata_file,
                  output_file,
                  class_name,
                  subclass_name=None,
                  where=None):
    """Create a text file which can be used as input for the LEfSe tool.

    This command
    1) collapses the input feature table at the genus level,
    2) computes relative frequency of the features,
    3) performs sample filtration if requested,
    4) changes the format of feature names,
    5) adds the relevant metadata as 'Class' and 'Subclass', and
    6) writes a text file which can be used as input for LEfSe.

    Parameters
    ----------
    table_file : str
        Path to the table file with the 'FeatureTable[Frequency]' type.
    taxonomy_file : str
        Path to the taxonomy file with the 'FeatureData[Taxonomy]' type.
    metadata_file : str
        Path to the metadata file.
    output_file : str
        Path to the output file.
    class_name : str
        Metadata column used as 'Class' by LEfSe.
    subclass_name : str, optional
        Metadata column used as 'Subclass' by LEfSe.
    where : str, optional
        SQLite 'WHERE' clause specifying sample metadata criteria.
    """
    _ = taxa.methods.collapse(
        table=Artifact.load(table_file),
        taxonomy=Artifact.load(taxonomy_file),
        level=6)

    _ = feature_table.methods.relative_frequency(
        table=_.collapsed_table)

    if where is None:
        df = _.relative_frequency_table.view(pd.DataFrame)
    else:
        _ = feature_table.methods.filter_samples(
            table=_.relative_frequency_table,
            metadata=Metadata.load(metadata_file),
            where=where)
        df = _.filtered_table.view(pd.DataFrame)

    def f(x):
        fields = x.split(';')
        for i in reversed(range(len(fields))):
            if fields[i] == '__':
                del fields[i]
            elif '__' in fields[i]:
                fields[i] = fields[i].split('__')[1]
            else:
                pass
        return '|'.join(fields)

    df.columns = [f(x) for x in df.columns.to_list()]

    mf = dokdo.get_mf(metadata_file)
    cols = mf.columns.to_list()
    df = pd.concat([df, mf], axis=1, join="inner")
    df.insert(0, "sample_id", df.index)
    df.insert(0, class_name, df.pop(class_name))
    cols.remove(class_name)

    if subclass_name is not None:
        df.insert(1, subclass_name, df.pop(subclass_name))
        cols.remove(subclass_name)

    df.drop(columns=cols, inplace=True)
    df.T.to_csv(output_file, header=False, sep='\t')
