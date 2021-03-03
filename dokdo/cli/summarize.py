import warnings
import pandas as pd
from qiime2 import Artifact

NO_VERBOSE_MESSAGE = ("There is no verbose option available "
                      "for the input Artifact file.")

def summarize(input_file, verbose=False):
    """Extract summary or verbose data from an Artifact file.

    This command automatically detects the input file's semantic type and
    then extracts summary or verbose data from it.

    Currently, the command supports the following semantic types:
        FeatureTable[Frequency]
        FeatureTable[RelativeFrequency]
        FeatureData[Sequence]
        FeatureData[AlignedSequence]

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
    elif str(artifact.type) in ["FeatureData[Sequence]",
        "FeatureData[AlignedSequence]"]:
        _parse_feature_data(artifact, verbose)
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

def _parse_feature_data(artifact, verbose):
    s = artifact.view(pd.Series)
    print("Number of features:", s.size)
    if verbose:
        print("Displaying only the first five records...")
        for index, value in s.apply(str).head().items():
            print(index)
            print(value)
