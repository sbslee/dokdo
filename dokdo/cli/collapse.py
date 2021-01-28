import os
import pandas as pd
from qiime2 import Artifact
from qiime2.plugins import taxa

def collapse(table_file, taxonomy_file, output_dir):
    os.mkdir(output_dir)
    for i in range(1, 8):
        _ = taxa.methods.collapse(
            table=Artifact.load(table_file),
            taxonomy=Artifact.load(taxonomy_file),
            level=i)
        _.collapsed_table.view(pd.DataFrame).T.to_csv(
            f"{output_dir}/level-{i}.csv")
