import pandas as pd

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
