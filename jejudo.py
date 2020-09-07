import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
import sklearn
import copy

def transform(jjd1, method):
    jjd2 = copy.deepcopy(jjd1)
    df = jjd2.asv_table

    if method == 'z':
        df = (df - df.mean()) / df.std()

    elif method == 'p':
        df = df / df.sum()

    else:
        raise ValueError("Incorrect method detected")

    jjd2.asv_table = df

    return jjd2

def ordinate(jejudo, method, n_components=2):
    udo = Udo()
    df = jejudo.asv_table.T

    if method == 'TSNE':
        embedding = sklearn.manifold.TSNE(n_components)
        X = embedding.fit_transform(df)

    elif method == 'PCA':
        embedding = sklearn.decomposition.PCA(n_components)
        X = embedding.fit_transform(df)

    elif method == 'MDS':
        embedding = sklearn.manifold.MDS(n_components)
        X = embedding.fit_transform(df)

    udo.method = method
    udo.embedding = embedding
    udo.X = X

    return udo

def plot_ordination(jejudo, udo, feature=None, elev=10, azim=30, figsize=(7,7)):
    n_components = udo.X.shape[1]

    if n_components == 2:
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_xlabel('1D')
        ax.set_ylabel('2D')
    elif n_components == 3:
        fig = plt.figure(1, figsize=figsize)
        ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=elev, azim=azim)
        ax.set_xlabel('1D')
        ax.set_ylabel('2D')
        ax.set_zlabel('3D')
    else:
        raise ValueError("Incorrect dimension detected")

    ax.set_title(f"{udo.method}")

    if feature:
        a = jejudo.smp_table[feature].unique()
        color_map = dict(zip(a, [None for x in a]))
        color_map = dict(sorted(color_map.items()))

        for k, v in color_map.items():
            i = jejudo.smp_table[feature] == k

            if n_components == 2:
                ax.scatter(udo.X[i, 0], udo.X[i, 1], label=k, color=v, s=90)
            else:
                ax.scatter(udo.X[i, 0], udo.X[i, 1], udo.X[i, 2], label=k, color=v, s=90)

        ax.legend(title=feature, loc="center left", bbox_to_anchor=(1, 0.5))

    else:
        if n_components == 2:
            ax.scatter(udo.X[:, 0], udo.X[:, 1], s=90)
        else:
            ax.scatter(udo.X[:, 0], udo.X[:, 1], udo.X[:, 2], s=90)

def plot_importance(jejudo, udo, n_tax=10):
    n_components = udo.embedding.components_.shape[0]

    a = {'ASV': jejudo.asv_table.index,
         '1D': udo.embedding.components_[0],
         '2D': udo.embedding.components_[1]}

    if n_components == 3:
        a['3D'] = udo.embedding.components_[2]

    df = pd.DataFrame(a)
    df = df.set_index('ASV')

    if n_components == 2:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,5))
    else:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(10,5))

    pc1 = df.abs().sort_values('1D', ascending=False)
    pc2 = df.abs().sort_values('2D', ascending=False)
    if n_components == 3:
        pc3 = df.abs().sort_values('3D', ascending=False)

    n_tax = 10

    ax1.bar(pc1.index[:n_tax], pc1['1D'][:n_tax])
    ax2.bar(pc2.index[:n_tax], pc2['2D'][:n_tax])
    if n_components == 3:
        ax3.bar(pc3.index[:n_tax], pc3['3D'][:n_tax])

    ax1.set_title('1D')
    ax2.set_title('2D')
    if n_components == 3:
        ax3.set_title('3D')

    ax1.set_ylabel('Feature Importance')

    ax1.tick_params(axis='x', rotation=90)
    ax2.tick_params(axis='x', rotation=90)
    if n_components == 3:
        ax3.tick_params(axis='x', rotation=90)

    fig.tight_layout()
    plt.show()

    ps1 = jejudo.tax_table.loc[ pc1.index[:n_tax] , : ]
    ps2 = jejudo.tax_table.loc[ pc2.index[:n_tax] , : ]
    if n_components == 3:
        ps3 = jejudo.tax_table.loc[ pc3.index[:n_tax] , : ]

    ps1.insert(0, 'Loading', df.loc[pc1.index[:n_tax], '1D'])
    ps2.insert(0, 'Loading', df.loc[pc2.index[:n_tax], '2D'])
    if n_components == 3:
        ps3.insert(0, 'Loading', df.loc[pc3.index[:n_tax], '3D'])

    display(ps1)
    display(ps2)
    if n_components == 3:
        display(ps3)

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

class Udo:
    def __init__(self):
        self._method = None
        self._embedding = None
        self._X = None

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, x):
        self._method = x

    @property
    def embedding(self):
        return self._embedding

    @embedding.setter
    def embedding(self, x):
        self._embedding = x

    @property
    def X(self):
        return self._X

    @X.setter
    def X(self, x):
        self._X = x

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