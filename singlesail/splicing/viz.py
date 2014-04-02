import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from singlesail.splicing.utils import get_switchy_score_order
from singlesail.colors import dark2
sns.set_axes_style('nogrid', 'talk')

def lavalamp(psi, color=None, title='', ax=None):
    """Make a 'lavalamp' scatter plot of many spliciang events

    Useful for visualizing many splicing events at once.

    Parameters
    ----------
    psi : array
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
        fig, ax = plt.subplots(figsize=(16,4))
    else:
        fig = plt.gcf()
    nrow, ncol = psi.shape
    x = np.vstack(np.ones(nrow) * i for i in range(ncol))

    color = dark2[0] if color is None else color

    try:
        # This is a pandas Dataframe
        # y = psi.values.T
        y = psi.values
    except AttributeError:
        # This is a numpy array
        # y = psi.T
        pass

    order = get_switchy_score_order(y)
    y = y[order, :]

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

    # Return the figure for saving
    return fig

def get_subplots_rows_cols(n):
    """Given a number of items that needs to be plotted, find the minimum
    number of rows and columns for subplot axes

    E.g. for plotting a number of cells' expression matrices.
    Totes stolen from pandas.Dataframe.hist() function.

    Parameters
    ----------
    n : int
        number of items to be plotted

    Returns
    -------
    nrows : int
        Number of rows for a matplotlib.pyplot.subplots() call
    ncols : int
        Number of columns for a matplotlib.pyplot.subplots() call
    """
    rows, cols = 1, 1
    while rows * cols < n:
        if cols > rows:
            rows += 1
        else:
            cols += 1
    return rows, cols