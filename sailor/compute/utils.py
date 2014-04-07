import numpy as np
import pandas as pd

def binify(df, binsize, vmin=0, vmax=1):
    """Makes a histogram of each row the provided binsize

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe whose rows you'd like to binify.
    binsize : float
        Size of bins
    vmin : float
        Minimum value of the bins
    vmax : float
        Maximum value of the bins

    Returns
    -------
    binned : pandas.DataFrame

    Raises
    ------


    """
    bins = np.arange(vmin, vmax + binsize, binsize)
    ncol = int(1.0 / binsize)
    nrow = df.shape[0]
    binned = np.zeros((nrow, ncol))

    # TODO: make sure this works for numpy matrices
    for i, (name, row) in enumerate(df.iterrows()):
        binned[i, :] = np.histogram(row, bins=bins, normed=True)[0]

    columns = ['{}-{}'.format(i, j) for i, j in zip(bins, bins[1:])]
    binned = pd.DataFrame(binned, index=df.index, columns=columns)
    return binned