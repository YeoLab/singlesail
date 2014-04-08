import numpy as np
from singlesail.data_model._Data import Data
from singlesail import compute
from sklearn.decomposition import PCA


class SplicingData(Data):
    def __init__(self, df, n_components=2, binsize=0.1,
                 reducer=PCA):
        """Instantiate a object for data scores with binned and reduced data

        Parameters
        ----------
        data : pandas.DataFrame
            A [n_events, n_samples] dataframe of splicing events
        n_components : int
            Number of components to use in the reducer
        binsize : float
            Value between 0 and 1, the bin size for binning the data scores
        reducer : sklearn.decomposition object
            An scikit-learn class that reduces the dimensionality of data
            somehow. Must accept the parameter n_components, have the
            functions fit, transform, and have the attribute components_

        """
        self.df = df
        self.reducer = reducer
        # self.psi_fillna_mean = self.data.T.fillna(self.data.mean(axis=1)).T
        self.binsize = binsize
        self.n_components = n_components
        self.binned = compute.utils.binify(self.df, binsize=self.binsize)
        self.reduce()

    def reduce(self):
        """Reduces dimensionality of the binned df score data
        """
        self.reducer_fit_to_binned = self.reducer(n_components=self.n_components).fit(self
                                                                    .binned)
        self.reduced_binned = self.reducer_fit_to_binned.transform(self.binned)
        if hasattr(self.reducer_fit_to_binned, 'explained_variance_ratio_'):
            self.plot_explained_variance(self.reducer_fit_to_binned,
                                         '{} on binned data'.format(self.reducer))
        return self
