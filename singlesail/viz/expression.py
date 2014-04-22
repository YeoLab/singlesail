__author__ = 'olga'
import matplotlib.pyplot as plt
import numpy as np

from singlesail.viz import add_subplot_axes

def coefficient_of_variation_plot(expression, log_base=2, ax=None, title='',
                                  plot_kws=None, ylim=None, xlim=None):
    """Plot the coefficient of variation

    Parameters
    ----------
    expression : pandas.DataFrame
        Expression count data, in the format (n_genes, n_samples). Assumed not
        log-normalized. If log-normalized, set log_base=1.
    log_base : float
        Which base log to use. Defaults to 2, for log2
    ax : matplotlib.axes.Axes
        Axes to plot on. Default None. If None, then grabs the current axes
    title :str
        Title of the

    Returns
    -------


    Raises
    ------
    """
    if ax is None:
        ax = plt.gca()

    expression = np.log(expression)/np.log(log_base)
    mean = expression.mean(axis=1)
    coefficient_of_variation = expression.std(axis=1) / mean

    if ylim is None:
        ymin = -1
        ymax = coefficient_of_variation.quantile(0.99)
    else:
        ymin, ymax = ylim

    if xlim is None:
        xmin = np.log(.01)


    plot_kws = {} if plot_kws is None else plot_kws
    plot_kws.setdefault('alpha', 0.25)
    plot_kws.setdefault('linewidth')
    plot_kws.setdefault('marker', 'o')

    ax.set_xscale('log', base=2)
    ax.plot(mean, coefficient_of_variation, **plot_kws)
    ax.set_title(title)
    ax.set_ylim(ymin, ymax)

    color = ax.lines[0]._get_markerfacecolor()

    subax = add_subplot_axes(ax, [0.7, 0.7, 0.3, 0.3])
    subax.hist(coefficient_of_variation,
               range=(0, coefficient_of_variation.max()), bins=50, color=color,
               log='y')
    subax.set_title('coeff of var')
    subax.set_xticklabels(map(int, subax.get_xticks()), rotation=60, )