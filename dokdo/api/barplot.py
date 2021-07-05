import tempfile
import dokdo
import pandas as pd
from .taxa_abundance_bar_plot import taxa_abundance_bar_plot
from qiime2 import Visualization
import matplotlib.pyplot as plt

def barplot(
    barplot_file, group, axis=0, figsize=(10, 10), level=1,
    count=0, items=None, by=None, label_columns=None,
    metadata=None, artist_kwargs=None, ylabel_fontsize=None,
    xaxis_repeated=False, cmap_name='Accent'
):
    """Create a grouped abundance bar plot.

    Under the hood, this method essentially wraps the
    `taxa_abundance_bar_plot` method.

    Parameters
    ----------
    barplot_file : str or qiime2.Visualization
        Visualization file or object from the q2-taxa plugin.
    group : str
        Metadata column.
    axis : int, default : 0
        By default, charts will be stacked vertically. Use 1 for horizontal
        stacking.
    figsize : tuple, default: (10, 10)
        Width, height in inches. Format: (float, float).
    level : int, default: 1
        Taxonomic level at which the features should be collapsed.
    count : int, default: 0
        The number of taxa to display. When 0, display all.
    items : list, optional
        Specify the order of charts.
    by : list, optional
        Column name(s) to be used for sorting the samples. Using 'index' will
        sort the samples by their name, in addition to other column name(s)
        that may have been provided. If multiple items are provided, sorting
        will occur by the order of the items.
    label_columns : list, optional
        The column names to be used as the x-axis labels.
    metadata : str or qiime2.Metadata, optional
        Metadata file or object.
    artist_kwargs : dict, optional
        Keyword arguments passed down to the _artist() method.
    ylabel_fontsize : float or str, optional
        Sets the y-axis label font size.
    xaxis_repeated : bool, default: False
        If true, remove all x-axis tick labels except for the bottom subplot.
        Ignored if `axis=1`.
    cmap_name : str, default: 'Accent'
        Name of the colormap passed to `matplotlib.cm.get_cmap()`.

    See Also
    --------
    dokdo.api.taxa_abundance_bar_plot

    Examples
    --------
    Below is a simple example.

    >>> barplot_file = f'{data_dir}/moving-pictures-tutorial/taxa-bar-plots.qzv'
    >>> dokdo.barplot(barplot_file,
    ...               'body-site',
    ...               axis=1,
    ...               figsize=(10, 6),
    ...               level=6,
    ...               count=8)

    .. image:: images/barplot-1.png

    We can draw the subplots vertically, which is particularly useful when the samples are matched.

    >>> dokdo.barplot(barplot_file,
    ...               'body-site',
    ...               axis=0,
    ...               figsize=(8, 10),
    ...               level=6,
    ...               count=8,
    ...               xaxis_repeated=True)

    .. image:: images/barplot-2.png
    """
    with tempfile.TemporaryDirectory() as t:
        vis = Visualization.load(barplot_file)
        vis.export_data(t)
        df = pd.read_csv(f'{t}/level-1.csv', index_col=0)

    if metadata is not None:
        mf = dokdo.get_mf(metadata)
        cols = _get_mf_cols(df)
        df.drop(columns=cols, inplace=True)
        df = pd.concat([df, mf], axis=1, join='inner')

    if items is None:
        _items = df[group].unique()
    else:
        _items = items

    if axis == 0:
        args = [len(_items), 3]
        gridspec_kw = dict(width_ratios=[0.01, 1, 0.01])
    else:
        args = [1, len(_items)+2]
        gridspec_kw=dict(width_ratios=[0.01]+[1 for x in _items]+[0.01])

    fig, axes = plt.subplots(*args, figsize=figsize, gridspec_kw=gridspec_kw)

    if artist_kwargs is None:
        artist_kwargs = {}

    _artist_kwargs = {'hide_ytexts': True, **artist_kwargs}

    plot_kwargs = dict(sort_by_mean2=False,
                       level=level,
                       count=count,
                       by=by,
                       label_columns=label_columns,
                       metadata=metadata,
                       cmap_name=cmap_name)

    if axis == 0:
        if xaxis_repeated:
            hide_xtexts = [True for x in range(len(axes[:, 1]))]
            hide_xtexts[-1] = False
        else:
            hide_xtexts = [False for x in range(len(axes[:, 1]))]

        for i, ax in enumerate(axes[:, 1]):
            taxa_abundance_bar_plot(barplot_file,
                                    ax=ax,
                                    include_samples={group: [_items[i]]},
                                    artist_kwargs={'title': _items[i],
                                                   'hide_xtexts': hide_xtexts[i],
                                                   **_artist_kwargs},
                                    **plot_kwargs)

    else:
        for i, ax in enumerate(axes[1:-1]):
            taxa_abundance_bar_plot(barplot_file,
                                    ax=ax,
                                    include_samples={group: [_items[i]]},
                                    artist_kwargs={'title': _items[i], **_artist_kwargs},
                                    **plot_kwargs)

    # Add the shared y-axis label.
    if axis == 0:
        gs = axes[0, 0].get_gridspec()
        for ax in axes[:, 0]:
            ax.remove()
        axbig = fig.add_subplot(gs[:, 0])
    else:
        axbig = axes[0]
    axbig.set_ylabel('Relative abundance (%)', fontsize=ylabel_fontsize)
    axbig.xaxis.set_visible(False)
    plt.setp(axbig.spines.values(), visible=False)
    axbig.tick_params(left=False, labelleft=False)
    axbig.patch.set_visible(False)

    # Add the shared legend.
    if axis == 0:
        gs = axes[0, -1].get_gridspec()
        for ax in axes[:, -1]:
            ax.remove()
        axbig = fig.add_subplot(gs[:, -1])
    else:
        axbig = axes[-1]

    taxa_abundance_bar_plot(barplot_file,
                            ax=axbig,
                            legend_short=True,
                            artist_kwargs={'legend_only': True,
                                           'legend_loc': 'center left',
                                           **_artist_kwargs},
                            **plot_kwargs)

    plt.tight_layout()
