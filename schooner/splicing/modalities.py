import numpy as np
import matplotlib.pyplot as plt


import brewer2mpl
from itertools import cycle
from scipy.spatial.distance import pdist, squareform
import skfuzzy
from sklearn.base import BaseEstimator, ClusterMixin, TransformerMixin
from sklearn.cluster import KMeans, spectral_clustering
from sklearn.decomposition import PCA
import seaborn as sns


class FuzzyCMeans(BaseEstimator, ClusterMixin, TransformerMixin):
    """Class for Fuzzy C-means clustering in scikit-learn cluster class format

    Implements the objective function as described:
    http://en.wikipedia.org/wiki/Fuzzy_clustering

    Parameters
    ----------
    n_clusters : int
        number of clusters to estimate
    exponent : float
        Exponent of objective function
    min_error : float
        The threshold at which if the objective function does not change by
        more than this, the iterations are stopped.
    max_iter : int
         Maximum number of iterations

    Attributes
    ----------
    `cluster_centers_` : array, [n_clusters, n_features]
        Coordinates of cluster centers

    `labels_` :
        Labels of each point

    `inertia_` : float
        The value of the inertia criterion associated with the chosen
        partition.

    Methods
    -------
    fit
        Like sklearn.cluster.KMeans, fit a matrix with samples to classify on
        the rows and features on the columns to the number of clusters.
        Basically, 'perform clustering'
    """
    def __init__(self, n_clusters=8, exponent=2, min_error=0.01, max_iter=100):
        """Initialize a Fuzzy C-Means clusterer

        Parameters
        ----------
        n_clusters : int
            number of clusters to estimate
        exponent : float
            Exponent of objective function
        min_error : float
            The threshold at which if the objective function does not change by
            more than this, the iterations are stopped.
        max_iter : int
             Maximum number of iterations

        Returns
        -------
        self : FuzzyCMeans
            An instantiated class of FuzzyCMeans, ready for clustering!
        """
        self.n_clusters = n_clusters
        self.exponent = exponent
        self.min_error = min_error
        self.max_iter = max_iter
        return self

    def fit(self, data, prob_thresh=None):
        """Fit the data to the number of clusters

        Parameters
        ----------
        data : numpy.array
            A numpy array of values (no NAs!) in the same format as required
            by scikit-learn, that in the shape [n_samples, n_features]

        Returns
        -------



        Raises
        ------

        """
        cluster_centers, fuzzy_matrix, initial_guess, distance_matrix, \
            objective_function_history, n_iter, fuzzy_partition_coeff = \
            skfuzzy.cmeans(data.T, c=self.n_clusters, m=self.exponent,
                          error=self.min_error, maxiter=self.max_iter)

        # rewrite as sklearn terminology
        self.cluster_centers_ = cluster_centers
        self.probability_of_labels = fuzzy_matrix
        self.labels_ = np.apply_along_axis(np.argmax, axis=0, arr=self.probability_of_labels)

        # Adjust labels for everything that didn't have a cluster membership prob >= 0.9
        self.unclassifiable = (self.probability_of_labels >= 0.9).sum(axis=0) == 0
        self.labels_[self.unclassifiable] = -1

        self.distance_matrix = distance_matrix
        self.objective_function_history = objective_function_history
        self.n_iter = n_iter
        self.fuzzy_partition_coeff = fuzzy_partition_coeff



class Data(object):
    def __init__(self, psi, n_components, step=0.1):
        self.psi = psi
        self.psi_fillna_mean = self.psi.T.fillna(self.psi.mean(axis=1)).T
        self.step = step
        self.n_components = n_components
        self.binify().reduce()

    def binify(self):
        self.bins = np.arange(0, 1+self.step, self.step)
        ncol = int(1/self.step)
        nrow = self.psi.shape[0]
        self.binned = np.zeros((nrow, ncol))
        for i, (name, row) in enumerate(self.psi.iterrows()):
            self.binned[i,:] = np.histogram(row, bins=self.bins, normed=True)[0]
        return self

    def reduce(self):
        self.pca_psi = PCA(n_components=self.n_components).fit(self.psi_fillna_mean)
        self.reduced_psi = self.pca_psi.transform(self.psi_fillna_mean)
        self.plot_explained_variance(self.pca_psi, 'PCA on psi (fillna with mean of event)')

        self.pca_binned = PCA(n_components=self.n_components).fit(self.binned)
        self.reduced_binned = self.pca_binned.transform(self.binned)
        self.plot_explained_variance(self.pca_binned, 'PCA on binned data')
        return self

    def plot_explained_variance(self, pca, title):
        # Plot the explained variance ratio
        fig, ax = plt.subplots()
        ax.plot(pca.explained_variance_ratio_, 'o-')
        ax.set_xticks(range(pca.n_components))
        ax.set_xticklabels(map(str, np.arange(pca.n_components)+1))
        ax.set_xlabel('Principal component')
        ax.set_ylabel('Fraction explained variance')
        ax.set_title(title)
        sns.despine()

    def calculate_distances(self, metric='euclidean'):
        self.pdist = squareform(pdist(self.binned, metric=metric))
        return self



def switchy_score(array):
    """Transform a 1D array of psi scores to a vector of "switchy scores"

    Calculates std deviation and mean of sine- and cosine-transformed
    versions of the array. Better than sorting by just the mean which doesn't
    push the really lowly variant events to the ends.

    Parameters
    ----------
    array : numpy.array
        A 1-D numpy array or something that could be cast as such (like a list)

    Returns
    -------
    float
        The "switchy score" of the data which can then be compared to other
        splicing event data

    @author Michael T. Lovci
    """
    array = np.array(array)
    variance = 1 - np.std(np.sin(array[~np.isnan(array)] * np.pi))
    mean_value = -np.mean(np.cos(array[~np.isnan(array)] * np.pi))
    return variance * mean_value

def get_switchy_score_order(x):
    """Apply switchy scores to a 2D array of psi scores

    Parameters
    ----------
    x : numpy.array
        A 2-D numpy array in the shape [n_events, n_samples]

    Returns
    -------
    numpy.array
        A 1-D array of the ordered indices, in switchy score order
    """
    switchy_scores = np.apply_along_axis(switchy_score, axis=0, arr=x)
    return np.argsort(switchy_scores)


class ClusteringTester(object):
    """Class for consistent evaluation of clustering methods

    Attributes
    ----------


    Methods
    -------
    hist_lavalamp
        Plot a histogram and lavalamp of psi scores from each cluster
    pca_viz
        Vizualize the clusters on the PCA of the data
    """
    def __init__(self, data, ClusterMethod, reduced='binned', cluster_kws=None,
                 colors=None):
        """Initialize ClusterTester and cluster the data

        Parameters
        ----------
        data : Data
            An object of the Data class
        ClusterMethod : sklearn.cluster class
            An object of the format from sklearn.cluster. Must have the fit()
            method, and create the attributes labels_
        reduced : str
            Specified which PCA-reduced data to use. Either the
            histogram-binned data ("binned") or the raw psi scores ("psi")
        """
        self.data = data
        self.reduced = self._get_reduced(reduced)

        cluster_kws = cluster_kws if cluster_kws is not None else {}
        self.clusterer = ClusterMethod(**cluster_kws)
        if ClusterMethod != spectral_clustering:
            self.clusterer.fit(self.reduced)
            self.labels = self.clusterer.labels_
        else:
            self.labels = self.clusterer
        self.labels_unique = set(self.labels)
        self.n_clusters = len(self.labels_unique)
        print 'n_clusters:', self.n_clusters

        if colors is None:
            self.colors = brewer2mpl.get_map('Set1', 'Qualitative', 8).mpl_colors
        else:
            self.colors = colors
        self.color_cycle = cycle(self.colors)

    def _get_reduced(self, reduced):
        """Sanely extracts the PCA-reduced data from the Data object

        Parameters
        ----------
        reduced : str
            Either "binned" or "psi"

        Returns
        -------
        reduced_data : numpy.array
            The PCA-reduced data from the specified array
        """
        if reduced.lower() == 'psi':
            reduced_data = self.data.reduced_psi
        elif reduced.lower() == 'binned':
            reduced_data = self.data.reduced_binned
        else:
            raise ValueError('Reduced data must be specified as one of "psi" '
                             'or "binned", not {}'.format(reduced))
        return reduced_data

    def _hist(self, ax, label, color):
        """Plot histograms of the psi scores of one label"""
        ax._hist(self.data.psi.ix[self.data.psi.index[self.labels == label],:].values.flat,
                bins=np.arange(0, 1.05, 0.05), facecolor=color, linewidth=0.1)
        ax.set_title('Cluster: {}'.format(label))
        ax.set_xlim(0,1)
        sns.despine()

    def _lavalamp(self, ax, label, color):
        """makes a _lavalamp of psi scores of one label"""
        nrow = (self.labels == label).sum()
        ncol = self.data.psi.shape[1]
        x = np.vstack(np.arange(nrow) for _ in range(ncol))
        y = self.data.psi.ix[self.data.psi.index[self.labels == label],:].values.T
        order = get_switchy_score_order(y)
        y = y[:,order]
        x_prev = x.max() + 1
        ax.scatter(x, y, color=color, alpha=0.5, edgecolor='#262626', linewidth=0.1)
        sns.despine()
        ax.set_xlim(0, x_prev)
        ax.set_ylim(0, 1)
        ax.set_title('n = {}'.format(nrow))

    def _annotate_centers(self, ax):
        """If the clusterer has cluster_centers_, plot the centroids

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
            Axes object to plot the annotation on
        """
        if type(self.clusterer) == KMeans or \
                type(self.clusterer == FuzzyCMeans):
            # Plot the centroids as a white X
            centroids = self.clusterer.cluster_centers_
        else:
            return
#             centroids = self.data.binned[self.clusterer.core_sample_indices_,:]

        ax.scatter(centroids[:, 0], centroids[:, 1],
                       marker='x', s=169, linewidths=3,
                       color='k', zorder=10)
        for i in range(self.n_clusters):
            try:
                ax.annotate(str(i),
                            (centroids[i,0], centroids[i,1]),
                            fontsize=24, xytext=(-20,0), textcoords='offset points')
            except:
                pass

    def hist_lavalamp(self):
        """Plot a histogram and _lavalamp for all the clusters

        Returns
        -------
        fig : matplotlib.pyplot.figure
            A figure instance with all the histograms and _lavalamp plots of
            all labels, for saving.
        """
        # Reset the color cycle in case we already cycled through it
        self.color_cycle = cycle(self.colors)

        fig = plt.figure(figsize=(16, 4*self.n_clusters))
        for i, (label, color) in enumerate(zip(self.labels_unique,
                                          self.color_cycle)):
            if label % 10 == 0:
                print 'plotting cluster {} of {}'.format(label, self.n_clusters)
            if label == -1:
                color = 'k'
            n_samples_in_cluster = (self.labels == label).sum()
            if n_samples_in_cluster <= 5:
                continue

            # fig = plt.figure(figsize=(16, 4))
            hist_ax = plt.subplot2grid((1, 5), (i,0), colspan=1, rowspan=1)
            lavalamp_ax = plt.subplot2grid((1,5), (i, 1), colspan=4, rowspan=4)

            self._hist(hist_ax, label, color=color)
            self._lavalamp(lavalamp_ax, label, color=color)
        return fig

    def pca_viz(self):
        """Visualizes the clusters on the PCA of the data

        Returns
        -------
        fig : matplotlib.pyplot.figure
            A figure instance with the PCA, for saving.

        """

        # Plot the decision boundary. For that, we will assign a color to each
        x_min, x_max = self.reduced[:, 0].min(), self.reduced[:, 0].max()
        y_min, y_max = self.reduced[:, 1].min(), self.reduced[:, 1].max()

        fig, ax = plt.subplots(figsize=(12,8))

        # Reset the color cycle in case we already cycled through it
        self.color_cycle = cycle(self.colors)

        colors = [color for _, color in zip(range(self.n_clusters), self.color_cycle)]
        color_list = [colors[int(label)] if label>=0 else 'k' for label in self.labels]
        ax.scatter(self.reduced[:, 0], self.reduced[:, 1],
                   color=color_list, alpha=0.25, linewidth=0.1, edgecolor='#262626')
        self._annotate_centers(ax)

        ax.set_title('{} clustering on the Motor Neuron dataset (PCA-reduced data)\n'
                 'Centroids are marked with black cross (step={:.2f})'.format(type(self.clusterer),
                                                                              self.data.step))
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_xticks(())
        ax.set_yticks(())
        sns.despine(left=True, bottom=True)
        return fig

    def violinplot_random_cluster_members(self, n=20):
        """Make violin plot of n random cluster members.

        Useful for seeing whether a cluster has bimodal events or not,
        which is not obvious from the lava lamp plot

        Parameters
        ----------
        n : int
            Number of cluster members to plot.

        Returns
        -------
        fig : matplotlib.pyplot.figure
            A figure instance with all the violin plots of all labels,
            for saving.

        """
        self.color_cycle = cycle(self.colors)
        fig, axes_array = plt.subplots(ncols=n,
                                       figsize=(n, 2*self.n_clusters),
                                       sharey=True)
        for axes, label, color in zip(axes_array, self.labels_unique,
                                self.color_cycle):
            if label == -1:
                color = 'k'
            these_labels = self.labels == label
            events = np.random.choice(self.data.psi.index[these_labels], size=n)
            y = self.data.psi.ix[events,:].values.T
            order = get_switchy_score_order(y)
            events = events[order]
            for event, ax in zip(events, axes):
        #         if i % 20 == 0:
                sns.violinplot(self.data.psi.ix[event], bw=0.1, inner='points',
                               color=color, linewidth=0, ax=ax, alpha=0.75)
                ax.set_ylim(0,1)
                ax.set_xticks([])
                ax.set_xlabel(label)
                sns.despine()
        return fig