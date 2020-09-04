import pandas as pd
import matplotlib.pyplot as plt
from sklearn import decomposition
from mpl_toolkits.mplot3d import Axes3D

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

    def pca2d(self, var):
        df = self.asv_table
        df = (df - df.mean()) / df.std()
        df = df.T

        fig, ax = plt.subplots(figsize=(7,7))
        pca = decomposition.PCA(2)
        X = pca.fit_transform(df)
        print(pca.explained_variance_ratio_)
        color_map = dict(zip(self.sample_data[var], self.sample_data.Color))

        for k, v in color_map.items():
            i = df.index.str.contains(k)
            ax.scatter(X[i, 0], X[i, 1], label=k, color=v, s=90)

        ax.set_xlabel("PC 1")
        ax.set_ylabel("PC 2")
        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.show()

    def pca3d(self, var, elev=10, azim=30):
        df = self.asv_table
        df = (df - df.mean()) / df.std()
        df = df.T

        fig = plt.figure(1, figsize=(7, 7))
        plt.clf()
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=elev, azim=azim)
        plt.cla()

        pca = decomposition.PCA(n_components=3)
        pca.fit(df)
        X = pca.transform(df)
        print(pca.explained_variance_ratio_)
        color_map = dict(zip(self.sample_data[var], self.sample_data.Color))

        for site, color in color_map.items():
            i = df.index.str.contains(site)
            ax.scatter(X[i, 0], X[i, 1], X[i, 2], label=site, color=color, s=80)

        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        ax.set_zlabel('PC 3')
        ax.w_xaxis.set_ticklabels([])
        ax.w_yaxis.set_ticklabels([])
        ax.w_zaxis.set_ticklabels([])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.show()