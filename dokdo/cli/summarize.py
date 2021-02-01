import pandas as pd
from qiime2 import Artifact

def summarize(input_file, verbose=False):
    """Extract summary or verbose data from an Artifact file.

    This command automatically detects the input file's semantic type and
    then extracts summary or verbose data from it.

    Parameters
    ----------
    input_file : str
        Path to the input Artifact file.
    verbose : bool, default: False
        Print a verbose version of the results.
    """
    artifact = Artifact.load(input_file)

    if str(artifact.type) in ["FeatureTable[Frequency]",
        "FeatureTable[RelativeFrequency]"]:
        _parse_feature_table(artifact, verbose)
    else:
        raise TypeError(f"Unsupported Artifact type: '{artifact.type}'")

def _parse_feature_table(artifact, verbose):
    df = artifact.view(pd.DataFrame)
    quantiles = [0, 0.25, 0.5, 0.75, 1]
    print("Number of samples:", df.shape[0])
    print("Number of features:", df.shape[1])
    print("Total frequency:", df.values.sum())
    print("Frequency per sample:")
    print(df.sum(axis=1).quantile(quantiles).to_string())
    print("Frequency per feature:")
    print(df.sum(axis=0).quantile(quantiles).to_string())
    if verbose:
        print("Samples:")
        print(" ".join(df.index.to_list()))
        print("Features:")
        print(" ".join(df.columns))
