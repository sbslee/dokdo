from qiime2 import Artifact

def summarize(input, verbose):
    table = Artifact.load(input)
    df = table.view(pd.DataFrame)
    quantiles = [0, 0.25, 0.5, 0.75, 1]
    print('Number of samples:', df.shape[0])
    print('Number of features:', df.shape[1])
    print('Total frequency:', df.values.sum())
    print('Frequency per sample:')
    print(df.sum(axis=1).quantile(quantiles).to_string())
    print('Frequency per feature:')
    print(df.sum(axis=0).quantile(quantiles).to_string())
    if verbose:
        print('Samples:')
        print(' '.join(df.index.to_list()))
        print('Features:')
        print(' '.join(df.columns))
