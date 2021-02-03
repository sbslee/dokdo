import os
import pandas as pd
from qiime2 import Artifact
from qiime2.plugins import taxa

def collapse(table_file, taxonomy_file, output_dir):
    """Create seven collapsed feature tables, one for each taxonomic
    level (i.e. 'level-1.csv' to 'level-7.csv').

    Parameters
    ----------
    table_file : str
        Path to the table file with the 'FeatureTable[Frequency]' type.
    taxonomy_file : str
        Path to the taxonomy file with the 'FeatureData[Taxonomy]' type.
    output_dir : str
        Path to the output directory.
    """
    os.mkdir(output_dir)
    for i in range(1, 8):
        _ = taxa.methods.collapse(
            table=Artifact.load(table_file),
            taxonomy=Artifact.load(taxonomy_file),
            level=i)
        _.collapsed_table.view(pd.DataFrame).T.to_csv(
            f"{output_dir}/level-{i}.csv")
