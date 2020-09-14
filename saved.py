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
