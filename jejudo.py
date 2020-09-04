import pandas as pd

def concat(jjd1, jjd2):
    jjd3 = Jejudo()

    df1 = jjd1.asv_table
    df2 = jjd2.asv_table
    df1["Sequence"] = jjd1.tax_table["Sequence"]
    df2["Sequence"] = jjd2.tax_table["Sequence"]
    df3 = pd.concat([df1, df2], ignore_index=True).fillna(0)
    asv_table = df3.groupby("Sequence").sum()

    df1 = jjd1.tax_table
    df2 = jjd2.tax_table
    df3 = pd.concat([df1, df2], ignore_index=True)
    df3 = df3.drop_duplicates(subset=["Sequence"])
    df3 = df3.sort_values("Sequence")
    tax_table = df3.set_index("Sequence")

    df1 = jjd1.sample_data
    df2 = jjd2.sample_data
    df1 = df1.reset_index()
    df2 = df2.reset_index()
    df3 = pd.concat([df1, df2], ignore_index=True)
    sample_data = df3.set_index("Sample")

    asv_list = [f"ASV{x+1}" for x in range(asv_table.shape[0])]
    asv_table.index = asv_list
    tax_table["Sequence"] = tax_table.index
    tax_table.index = asv_list

    jjd3.asv_table = asv_table
    jjd3.tax_table = tax_table
    jjd3.sample_data = sample_data

    return jjd3

class Jejudo:
    def __init__(self):
        self._asv_table = pd.DataFrame()
        self._tax_table = pd.DataFrame()
        self._sample_data = pd.DataFrame()

    @property
    def asv_table(self):
        return self._asv_table

    @asv_table.setter
    def asv_table(self, x):
        self._asv_table = x

    @property
    def tax_table(self):
        return self._tax_table

    @tax_table.setter
    def tax_table(self, x):
        self._tax_table = x

    @property
    def sample_data(self):
        return self._sample_data

    @sample_data.setter
    def sample_data(self, x):
        self._sample_data = x

    def describe(self):
        print("asv_table", "ASV Table:", self._asv_table.shape)
        print("tax_table", "Taxonomy Table:", self._tax_table.shape)
        print("sample_data", "Sample Data:", self._sample_data.shape)

    def read_files(self, asv_table_file, tax_table_file, sample_data_file):
        asv_table = pd.read_csv(asv_table_file, index_col=0)
        tax_table = pd.read_csv(tax_table_file, index_col=0)
        sample_data = pd.read_table(sample_data_file, index_col=0)

        asv_table = asv_table.T
        asv_list = [f"ASV{x+1}" for x in range(asv_table.shape[0])]
        asv_table.index = asv_list
        tax_table["Sequence"] = tax_table.index
        tax_table.index = asv_list

        self.asv_table = asv_table
        self.tax_table = tax_table
        self.sample_data = sample_data

    def collapse(self, rank):
        ranks = ["Kingdom", "Phylum", "Class", "Order",
                 "Family", "Genus", "Species"]
        i = ranks.index(rank) + 1
        df = self.asv_table
        a = self.tax_table.iloc[:, :i].astype(str)
        df["Target"] = a.agg(":".join, axis=1)
        df = df.groupby("Target").sum()
        return df
