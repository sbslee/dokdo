import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans
import sklearn
import copy
from matplotlib_venn import venn2
import matplotlib
import seaborn as sns

TAXA = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

def shannon_diversity_index(l):
    def p(n, N):
        return (n/N) * np.log(n/N)
    N = sum(l)
    return -sum(p(n, N) for n in l if n)

def simpson_diversity_index(l):
    def p(n, N):
        return (n/N)**2
    N = sum(l)
    return sum(p(n, N) for n in l if n)

def inverse_simpson_diversity_index(l):
    return 1 / simpson_diversity_index(l)

def plot_richness(jjd, feature, measure='Shannon', ax=None, show=False,
                  figsize=None):
    measures = {'Shannon': shannon_diversity_index,
                'Simpson': shannon_diversity_index,
                'InverseSimpson': inverse_simpson_diversity_index}

    if not ax:
        fig, ax = plt.subplots(figsize=figsize)

    df = jjd.asv_table
    a = jjd.smp_table[feature]
    b = df.apply(measures[measure], axis=0)
    df = pd.DataFrame({feature: a, measure: b})
    sns.boxplot(x=feature, y=measure, data=df, ax=ax)
    ax.set_ylabel(f"Alpha Diversity Measure ({measure})")

    if show:
        plt.show()

def plot_summary(jjd, method, feature=None, ax=None, show=False,
                 figsize=(7,7), bins=None, group='Bacteria'):
    if not ax:
        fig, ax = plt.subplots(figsize=figsize)

    if method == 'fh': # Feature histogram
        ax.hist(jjd.smp_table[feature])
        ax.set_xlabel(feature)
        ax.set_ylabel('Count')

    elif method == 'nz': # ASVs with non-zero values
        a = jjd.asv_table.astype(bool).sum(axis=1)
        ax.hist(a, bins=bins)
        ax.set_xlabel("Number of Non-Zero Values")
        ax.set_ylabel("Number of ASVs")

    elif method == 'rl': # Read length
        a = jjd.seq_table['Sequence'].str.len()
        ax.hist(a, bins=bins)
        ax.set_xlabel("Read Length (bp)")
        ax.set_ylabel("Number of ASVs")

    elif method == 'nr': # Number of reads
        d = {feature: jjd.smp_table[feature],
             'ASV': jjd.asv_table.T.sum(axis=1) / 1e3}
        df = pd.DataFrame(d)
        sns.boxplot(x=feature, y="ASV", data=df, ax=ax)
        ax.set_xlabel(feature)
        ax.set_ylabel("Number of Reads (K)")

    elif method == 'nt': # Number of taxa by group
        df = collapse(jjd, 'Kingdom', method='u').T
        df[feature] = jjd.smp_table[feature]
        sns.boxplot(x=feature, y=group, data=df, ax=ax)
        ax.set_xlabel(feature)
        ax.set_ylabel(f"Number of taxa ({group})")

    else:
        raise ValueError("Incorrect method detected")

    if show:
        plt.show()

def collapse(jjd, rank, method='c'):
    if method == 'c':
        df = _collapse_counts(jjd, rank)
    elif method == 'u':
        df = _collapse_unique(jjd, rank)
    else:
        raise ValueError("Incorrect method detected")

    return df

def _collapse_counts(jjd, rank):
    i = TAXA.index(rank) + 1
    df = jjd.asv_table
    a = jjd.tax_table.iloc[:, :i].astype(str)
    df['Target'] = a.agg(':'.join, axis=1)
    df = df.groupby('Target').sum()
    return df

def _collapse_unique(jjd, rank):
    i = TAXA.index(rank) + 1
    df = jjd.asv_table.astype(bool)
    a = jjd.tax_table.iloc[:, :i].astype(str)
    df['Target'] = a.agg(':'.join, axis=1)
    df = df.groupby('Target').sum()
    return df

def subset(jjd1, var, target):
    i = jjd1.smp_table[var] == target
    smp_df = jjd1.smp_table.loc[i, :]
    asv_df = jjd1.asv_table.loc[:, i]

    i = asv_df.sum(axis=1) != 0
    asv_df = asv_df.loc[i, :]
    tax_df = jjd1.tax_table.loc[i, :]
    seq_df = jjd1.seq_table.loc[i, :]

    asv_list = [f"ASV{x+1}" for x in range(asv_df.shape[0])]
    asv_df.index = asv_list
    tax_df.index = asv_list
    seq_df.index = asv_list

    jjd2 = Jejudo()
    jjd2.asv_table = asv_df
    jjd2.tax_table = tax_df
    jjd2.smp_table = smp_df
    jjd2.seq_table = seq_df

    return jjd2

def remove(jjd1, method, n_samples=0, n_bases=1000):
    jjd2 = copy.deepcopy(jjd1)

    asv_df = jjd1.asv_table
    tax_df = jjd1.tax_table
    seq_df = jjd1.seq_table

    if method == 's':
        i = jjd1.asv_table.astype(bool).sum(axis=1) >= n_samples
    elif method == 'l':
        i = jjd1.seq_table['Sequence'].str.len() >= n_bases
    else:
        raise ValueError("Incorrect method detected")

    asv_df = asv_df[i]
    tax_df = tax_df[i]
    seq_df = seq_df[i]

    asv_list = [f"ASV{x+1}" for x in range(asv_df.shape[0])]
    asv_df.index = asv_list
    tax_df.index = asv_list
    seq_df.index = asv_list

    jjd2.asv_table = asv_df
    jjd2.tax_table = tax_df
    jjd2.seq_table = seq_df
    
    return jjd2

def plot_comparison(jjd1, jjd2, figsize=(12,6), fontsize=20):
    df1 = jjd1.asv_table
    df2 = jjd2.asv_table

    if not all(df1.columns == df2.columns):
        raise ValueError("Samples not matched")

    df1['Sequence'] = jjd1.seq_table['Sequence']
    df2['Sequence'] = jjd2.seq_table['Sequence']

    df3 = pd.merge(df1, df2, how='inner', on=['Sequence'])
    del df3['Sequence']

    n = int(df3.shape[1] / 2)

    df1 = df3.iloc[:, :n].stack()
    df2 = df3.iloc[:, n:].stack()

    r2 = sklearn.metrics.r2_score(df1, df2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    AB = df3.shape[0]
    Ab = jjd1.asv_table.shape[0] - AB
    aB = jjd2.asv_table.shape[0] - AB
    vd = venn2(subsets=(Ab, aB, AB), ax=ax1)

    for x in vd.set_labels:
        x.set_fontsize(fontsize)

    for x in vd.subset_labels:
        x.set_fontsize(fontsize)

    ax2.set_xlabel('A', fontsize=fontsize)
    ax2.set_ylabel('B', fontsize=fontsize)
    ax2.plot([0, 1], [0, 1], color='red', transform=ax2.transAxes)
    ax2.scatter(df1, df2)
    ax2.text(0.7, 0.1, 'R-squared = %0.4f' % r2, fontsize=fontsize,
            horizontalalignment='center', verticalalignment='center',
            transform = ax2.transAxes)

def transform(jjd1, method):
    jjd2 = copy.deepcopy(jjd1)
    df = jjd2.asv_table

    if method == 'z':
        df = (df - df.mean()) / df.std()
    elif method == 'p':
        df = df / df.sum()
    elif method == 's':
        df = df.pow(1/2)
    elif method == 'l':
        df = np.log(df + 0.1)
    else:
        raise ValueError("Incorrect method detected")

    jjd2.asv_table = df

    return jjd2

def ordinate(jjd, method, n_components=2):
    ud = Udo()
    df = jjd.asv_table.T

    if method == 'TSNE':
        embedding = sklearn.manifold.TSNE(n_components)
        X = embedding.fit_transform(df)

    elif method == 'PCA':
        embedding = sklearn.decomposition.PCA(n_components)
        X = embedding.fit_transform(df)

    elif method == 'MDS':
        embedding = sklearn.manifold.MDS(n_components)
        X = embedding.fit_transform(df)

    elif method == 'NMDS':
        embedding = sklearn.manifold.MDS(n_components, metric=False)
        X = embedding.fit_transform(df)

    else:
        raise ValueError("Incorrect method detected")

    ud.method = method
    ud.embedding = embedding
    ud.X = X

    return ud

def plot_ordination(jjd, ud, *args, **kwargs):

    n_components = ud.X.shape[1]

    if n_components == 2:
        _plot_2d_ordination(jjd, ud, *args, **kwargs)
    elif n_components == 3:
        _plot_3d_ordination(jjd, ud, *args, **kwargs)
    else:
        raise ValueError("Incorrect dimension detected")

def _plot_2d_ordination(jjd, ud, ax=None, show=False, s=20, color=None, 
                        figsize=None, legend=True, labels=True, ticks=True, 
                        title=True):
    if not ax:
        fig, ax = plt.subplots(figsize=figsize)

    if title:
        ax.set_title(f"{ud.method}")

    if not ticks:
        ax.set_xticks([])
        ax.set_yticks([])

    if labels:
        ax.set_xlabel('1D')
        ax.set_ylabel('2D')

    if color:
        for c in sorted(jjd.smp_table[color].unique()):
            i = jjd.smp_table[color] == c
            ax.scatter(ud.X[i, 0], ud.X[i, 1], s=s, label=c)

        if legend:
            ax.legend(title=color, loc="center left", bbox_to_anchor=(1, 0.5))

    else:
        ax.scatter(ud.X[:, 0], ud.X[:, 1])

    if show:
        plt.show()

def _plot_3d_ordination(jjd, ud, ax=None, show=False, s=20, color=None,
                        figsize=None, legend=True, labels=True, ticks=True,
                        title=True, elev=30, azim=30):
    if not ax:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.view_init(elev=elev, azim=azim)

    if title:
        ax.set_title(f"{ud.method}")

    if not ticks:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])

    if labels:
        ax.set_xlabel('1D')
        ax.set_ylabel('2D')
        ax.set_zlabel('3D')

    if color:
        for c in sorted(jjd.smp_table[color].unique()):
            i = jjd.smp_table[color] == c
            ax.scatter(ud.X[i, 0], ud.X[i, 1], ud.X[i, 2], s=s, label=c)

        if legend:
            ax.legend(title=color, loc="center left", bbox_to_anchor=(1, 0.5))

    else:
        ax.scatter(ud.X[:, 0], ud.X[:, 1], ud.X[:, 2])

    if show:
        plt.show()

def plot_importance(jjd, ud, n_tax=10):
    n_components = ud.embedding.components_.shape[0]

    a = {'ASV': jjd.asv_table.index,
         '1D': ud.embedding.components_[0],
         '2D': ud.embedding.components_[1]}

    if n_components == 3:
        a['3D'] = ud.embedding.components_[2]

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

    ps1 = jjd.tax_table.loc[ pc1.index[:n_tax] , : ]
    ps2 = jjd.tax_table.loc[ pc2.index[:n_tax] , : ]
    if n_components == 3:
        ps3 = jjd.tax_table.loc[ pc3.index[:n_tax] , : ]

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

    asv_df = asv_df.sort_index(axis=1)
    smp_df = smp_df.sort_index()

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

    def kmeans(self, n):
        df = self.asv_table
        df = (df - df.mean()) / df.std()
        X = df.T
        kmeans = KMeans(n_clusters=n)
        kmeans.fit(X)
        y = kmeans.predict(X)
        self.smp_table = self.smp_table.assign(KMeans=y)