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
                  class_col,
                  subclass_col=None,
                  subject_col=None,
                  where=None):
    """Create a TSV file which can be used as input for the LEfSe tool.

    This command
    1) collapses the input feature table at the genus level,
    2) computes relative frequency of the features,
    3) performs sample filtration if requested,
    4) changes the format of feature names,
    5) adds the relevant metadata as 'Class', 'Subclass', and 'Subject', and
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
    class_col : str
        Metadata column used as 'Class' by LEfSe.
    subclass_col : str, optional
        Metadata column used as 'Subclass' by LEfSe.
    subject_col : str, optional
        Metadata column used as 'Subject' by LEfSe.
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
        for c in ['-', '[', ']', '(', ')', ' ']:
            x = x.replace(c, '_')

        ranks = x.split(';')
        base = ranks[0]
        result = [base]

        for i, rank in enumerate(ranks[1:], start=2):
            if rank == '__':
                result.append(f'{base}_x__L{i}')
            elif rank.split('__')[1] == '':
                result.append(f'{base}_{rank}L{i}')
            else:
                result.append(rank)
                base = rank

        return '|'.join(result)

    df.columns = [f(x) for x in df.columns.to_list()]

    mf = dokdo.get_mf(metadata_file)
    mf = mf.replace(' ', '_', regex=True)
    cols = mf.columns.to_list()
    df = pd.concat([df, mf], axis=1, join="inner")
    df.insert(0, class_col, df.pop(class_col))
    cols.remove(class_col)

    if subclass_col is None and subject_col is None:
        pass
    elif subclass_col is not None and subject_col is None:
        df.insert(1, subclass_col, df.pop(subclass_col))
        cols.remove(subclass_col)
    elif subclass_col is None and subject_col is not None:
        df.insert(1, subject_col, df.pop(subject_col))
        cols.remove(subject_col)
    else:
        df.insert(1, subclass_col, df.pop(subclass_col))
        df.insert(2, subject_col, df.pop(subject_col))
        cols.remove(subclass_col)
        cols.remove(subject_col)

    df.drop(columns=cols, inplace=True)
    df.T.to_csv(output_file, header=False, sep='\t')
