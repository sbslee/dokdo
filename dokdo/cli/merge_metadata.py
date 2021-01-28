def merge_metadata(metadata,
                   output):
    dfs = []

    for file in metadata:
        dfs.append(Metadata.load(file).to_dataframe())

    Metadata(pd.concat(dfs)).save(output)
