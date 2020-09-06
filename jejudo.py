import pandas as pd
import matplotlib.pyplot as plt
from sklearn import decomposition
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans

def concat(jjd1, jjd2):
    jjd3 = Jejudo()

    df1 = jjd1.asv_table
    df2 = jjd2.asv_table
    df1['Sequence'] = jjd1.seq_table['Sequence']
    df2['Sequence'] = jjd2.seq_table['Sequence']
    df3 = pd.concat([df1, df2], ignore_index=True).fillna(0)
    df3 = df3.groupby('Sequence').sum()
    asv_df = df3.sort_values('Sequence')

    df1 = jjd1.tax_table
    df2 = jjd2.tax_table
    df1['Sequence'] = jjd1.seq_table['Sequence']
    df2['Sequence'] = jjd2.seq_table['Sequence']
    df3 = pd.concat([df1, df2], ignore_index=True)
    df3 = df3.drop_duplicates(subset=['Sequence'])
    df3 = df3.sort_values('Sequence')
    tax_df = df3.set_index('Sequence')

    if sum(asv_df.index == tax_df.index) != asv_df.shape[0]:
        raise ValueError("Sequences not matched")

    df1 = jjd1.smp_table
    df2 = jjd2.smp_table
    df1 = df1.reset_index()
    df2 = df2.reset_index()
    df3 = pd.concat([df1, df2], ignore_index=True)
    smp_df = df3.set_index('Sample')

    asv_list = [f"ASV{x+1}" for x in range(asv_df.shape[0])]
    asv_df.index = asv_list
    d = {'ASV': asv_list, 'Sequence': tax_df.index}
    seq_df = pd.DataFrame(d)
    seq_df = seq_df.set_index('ASV')
    tax_df.index = asv_list

    jjd3.asv_table = asv_df
    jjd3.tax_table = tax_df
    jjd3.smp_table = smp_df
    jjd3.seq_table = seq_df

    return jjd3

class Jejudo:
    def __init__(self):
        self._asv_table = pd.DataFrame()
        self._tax_table = pd.DataFrame()
        self._smp_table = pd.DataFrame()
        self._seq_table = pd.DataFrame()

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
    def smp_table(self):
        return self._smp_table

    @smp_table.setter
    def smp_table(self, x):
        self._smp_table = x

    @property
    def seq_table(self):
        return self._seq_table

    @seq_table.setter
    def seq_table(self, x):
        self._seq_table = x

    def describe(self):
        a1, a2 = self._asv_table.shape
        b1, b2 = self._tax_table.shape
        c1, c2 = self._smp_table.shape
        d1, d2 = self._seq_table.shape

        print('asv_table:', f'[ {a1} taxa and {a2} samples ]')
        print('tax_table:', f'[ {b1} taxa and {b2} ranks ]')
        print('smp_table:', f'[ {c1} samples and {c2} features ]')
        print('seq_table:', f'[ {d1} sequences ]')

    def read_files(self, asv_file, tax_file, smp_file):
        asv_df = pd.read_csv(asv_file, index_col=0)
        tax_df = pd.read_csv(tax_file, index_col=0)
        smp_df = pd.read_table(smp_file, index_col=0)

        asv_df = asv_df.sort_index()
        smp_df = smp_df.sort_index()

        if not all(asv_df.index == smp_df.index):
            raise ValueError("Samples differ in ASV Table and Sample Table")

        asv_df = asv_df.T
        asv_list = [f"ASV{x+1}" for x in range(asv_df.shape[0])]
        asv_df.index = asv_list
        d = {'ASV': asv_list, 'Sequence': tax_df.index}
        seq_df = pd.DataFrame(d)
        seq_df = seq_df.set_index('ASV')
        tax_df.index = asv_list

        self.asv_table = asv_df
        self.tax_table = tax_df
        self.smp_table = smp_df
        self.seq_table = seq_df

    def collapse(self, rank):
        ranks = ['Kingdom', 'Phylum', 'Class', 'Order',
                 'Family', 'Genus', 'Species']
        i = ranks.index(rank) + 1
        df = self.asv_table
        a = self.tax_table.iloc[:, :i].astype(str)
        df['Target'] = a.agg(':'.join, axis=1)
        df = df.groupby('Target').sum()
        return df

    def pca2d(self, var, color_map=None):
        result = {}

        df = self.asv_table
        df = (df - df.mean()) / df.std()
        df = df.T

        fig, ax = plt.subplots(figsize=(7,7))
        pca = decomposition.PCA(2)
        X = pca.fit_transform(df)

        # Store the percentage of variance explained by each PC.
        a = pca.explained_variance_ratio_

        summary = ["Variance explained by PC 1: %.2f" % a[0],
                   "Variance explained by PC 2: %.2f" % a[1],
                   "Total variance explained: %.2f" % (a[0] + a[1])]

        result['summary'] = '\n'.join(summary)

        # Store the PCA plot.
        if not color_map:
            a = self.smp_table[var].unique()
            color_map = dict(zip(a, [None for x in a]))

        color_map = dict(sorted(color_map.items()))

        for k, v in color_map.items():
            i = self.smp_table[var] == k
            ax.scatter(X[i, 0], X[i, 1], label=k, color=v, s=90)

        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.close()
        result['plot1'] = fig

        # Store the feature importance plot.
        fig, (ax1, ax2) = plt.subplots(1, 2)

        a = {'ASV': df.columns, 'PC1': pca.components_[0], 'PC2': pca.components_[1]}
        df = pd.DataFrame(a)
        df = df.set_index('ASV')

        pc1 = df.abs().sort_values('PC1', ascending=False)
        pc2 = df.abs().sort_values('PC2', ascending=False)

        n = 10

        ax1.bar(pc1.index[:n], pc1['PC1'][:n])
        ax2.bar(pc2.index[:n], pc2['PC2'][:n])

        ax1.set_title('PC 1')
        ax2.set_title('PC 2')

        ax1.set_ylabel('Feature Importance')
        ax2.set_ylabel('Feature Importance')

        ax1.tick_params(axis='x', rotation=90)
        ax2.tick_params(axis='x', rotation=90)

        fig.tight_layout()
        plt.close()
        result['plot2'] = fig

        # Store the top ASV loadings for each PC.
        pd.set_option('display.max_colwidth', None)

        ps1 = self.tax_table.loc[ pc1.index[:n] , : ]
        ps2 = self.tax_table.loc[ pc2.index[:n] , : ]

        ps1.insert(0, 'Loading', df.loc[pc1.index[:n], 'PC1'])
        ps2.insert(0, 'Loading', df.loc[pc2.index[:n], 'PC2'])

        result['loading1'] = ps1
        result['loading2'] = ps2

        return result

    def pca3d(self, var, color_map=None, elev=10, azim=30):
        result = {}

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

        # Store the percentage of variance explained by each PC.
        a = pca.explained_variance_ratio_

        summary = ["Variance explained by PC 1: %.2f" % a[0],
                   "Variance explained by PC 2: %.2f" % a[1],
                   "Variance explained by PC 3: %.2f" % a[2],
                   "Total variance explained: %.2f" % (a[0] + a[1] + a[2])]

        result['summary'] = '\n'.join(summary)

        # Store the PCA plot.
        if not color_map:
            a = self.smp_table[var].unique()
            color_map = dict(zip(a, [None for x in a]))

        color_map = dict(sorted(color_map.items()))

        for k, v in color_map.items():
            i = self.smp_table[var] == k
            ax.scatter(X[i, 0], X[i, 1], X[i, 2], label=k, color=v, s=80)

        ax.set_xlabel('PC 1')
        ax.set_ylabel('PC 2')
        ax.set_zlabel('PC 3')
        ax.w_xaxis.set_ticklabels([])
        ax.w_yaxis.set_ticklabels([])
        ax.w_zaxis.set_ticklabels([])
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        plt.close()
        result['plot1'] = fig

        # Store the feature importance plot.
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)

        a = {'ASV': df.columns, 'PC1': pca.components_[0],
             'PC2': pca.components_[1], 'PC3': pca.components_[2]}
        df = pd.DataFrame(a)
        df = df.set_index('ASV')

        pc1 = df.abs().sort_values('PC1', ascending=False)
        pc2 = df.abs().sort_values('PC2', ascending=False)
        pc3 = df.abs().sort_values('PC3', ascending=False)

        n = 10

        ax1.bar(pc1.index[:n], pc1['PC1'][:n])
        ax2.bar(pc2.index[:n], pc2['PC2'][:n])
        ax3.bar(pc3.index[:n], pc3['PC3'][:n])

        ax1.set_title('PC 1')
        ax2.set_title('PC 2')
        ax3.set_title('PC 3')

        ax1.set_ylabel('Feature Importance')
        ax2.set_ylabel('Feature Importance')
        ax3.set_ylabel('Feature Importance')

        ax1.tick_params(axis='x', rotation=90)
        ax2.tick_params(axis='x', rotation=90)
        ax3.tick_params(axis='x', rotation=90)

        fig.tight_layout()
        plt.close()
        result['plot2'] = fig

        # Store the top ASV loadings for each PC.
        pd.set_option('display.max_colwidth', None)

        ps1 = self.tax_table.loc[ pc1.index[:n] , : ]
        ps2 = self.tax_table.loc[ pc2.index[:n] , : ]
        ps3 = self.tax_table.loc[ pc3.index[:n] , : ]

        ps1.insert(0, 'Loading', df.loc[pc1.index[:n], 'PC1'])
        ps2.insert(0, 'Loading', df.loc[pc2.index[:n], 'PC2'])
        ps3.insert(0, 'Loading', df.loc[pc3.index[:n], 'PC3'])

        result['loading1'] = ps1
        result['loading2'] = ps2
        result['loading3'] = ps3

        return result

    def subset(self, var, target):
        i = self.smp_table[var] == target
        smp_df = self.smp_table.loc[i, :]
        asv_df = self.asv_table.loc[:, i]

        i = asv_df.sum(axis=1) != 0
        asv_df = asv_df.loc[i, :]
        tax_df = self.tax_table.loc[i, :]
        seq_df = self.seq_table.loc[i, :]

        asv_list = [f"ASV{x+1}" for x in range(asv_df.shape[0])]
        asv_df.index = asv_list
        tax_df.index = asv_list
        seq_df.index = asv_list

        jjd = Jejudo()
        jjd.asv_table = asv_df
        jjd.tax_table = tax_df
        jjd.smp_table = smp_df
        jjd.seq_table = seq_df

        return jjd

    def kmeans(self, n):
        df = self.asv_table
        df = (df - df.mean()) / df.std()
        X = df.T
        kmeans = KMeans(n_clusters=n)
        kmeans.fit(X)
        y = kmeans.predict(X)
        self.smp_table = self.smp_table.assign(KMeans=y)