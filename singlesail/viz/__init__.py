__author__ = 'olga'

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from singlesail.colors import dark2
from singlesail.splicing.utils import get_switchy_score_order

sns.set(style='white', context='talk')


def lavalamp(psi, color=None, title='', ax=None):
    """Make a 'lavalamp' scatter plot of many spliciang events

    Useful for visualizing many splicing events at once.

    Parameters
    ----------
    df : array
        A (n_events, n_samples) matrix either as a numpy array or as a pandas
        DataFrame

    color : matplotlib color
        Color of the scatterplot. Defaults to a dark teal

    title : str
        Title of the plot. Default ''

    ax : matplotlib.Axes object
        The axes to plot on. If not provided, will be created


    Returns
    -------
    fig : matplotlib.Figure
        A figure object for saving.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 4))
    else:
        fig = plt.gcf()
    nrow, ncol = psi.shape
    x = np.vstack(np.arange(nrow) for _ in range(ncol))

    color = dark2[0] if color is None else color

    try:
        # This is a pandas Dataframe
        y = psi.values.T
    except AttributeError:
        # This is a numpy array
        y = psi.T

    order = get_switchy_score_order(y)
    y = y[:, order]

    # Add one so the last value is actually included instead of cut off
    xmax = x.max() + 1
    ax.scatter(x, y, color=color, alpha=0.5, edgecolor='#262626', linewidth=0.1)
    sns.despine()
    ax.set_ylabel('$\Psi$')
    ax.set_xlabel('{} splicing events'.format(nrow))
    ax.set_xticks([])

    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 1)
    ax.set_title(title)


def reducedplot(reduced, color=dark2[0],
                x_dim=0, y_dim=1, title='', ax=None,
                groupby=None, group_to_color=None):
    """Visualizes reduced data

    Returns
    -------
    fig : matplotlib.pyplot.figure
        A figure instance with the PCA, for saving.

    """
    X = reduced[:, x_dim]
    Y = reduced[:, y_dim]
    x_min, x_max = X.min(), X.max()
    y_min, y_max = Y.min(), Y.max()

    if ax is None:
        ax = plt.gca()

    if groupby is None:
        ax.scatter(X, Y,
                    color=color,
                   alpha=0.25, linewidth=0.1, edgecolor='#262626')
    else:
        for name, df in reduced.groupby(groupby):
            X = df[:, x_dim]
            Y = df[:, y_dim]
            color = group_to_color[name]
            ax.scatter(X, Y, color=color, alpha=0.25, linewidth=0.1,
                       edgecolor='#262626', label=name)
        ax.legend()

    ax.set_title(title)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xticks(())
    ax.set_yticks(())
    sns.despine(left=True, bottom=True)

def print_hello():
    print "hello"