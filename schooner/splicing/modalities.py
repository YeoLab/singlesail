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
        self.binify()
        self.reduce()
        return self

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
    """
    Transforms a psi score (between 0 and 1) and provides a reasonable way to
    sort by psi scores on lava lamp plots.

    array is something that can be cast to a 1-D np array

    @author Michael T. Lovci
    """
    array = np.array(array)
    variance = 1 - np.std(np.sin(array[~np.isnan(array)] * np.pi))
    mean_value = -np.mean(np.cos(array[~np.isnan(array)] * np.pi))
    return variance * mean_value

def get_switchy_score_order(x):
    switchy_scores = np.apply_along_axis(switchy_score, axis=0, arr=x)
    return np.argsort(switchy_scores)

class ClusteringTester(object):
    def __init__(self, data, ClusterMethod, reduced='binned', cluster_kws=None, colors=None):
        self.data = data
        self.reduced = self.data.reduced_psi if reduced is 'psi' else self.data.reduced_binned
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

    def hist(self, ax, label, color):
        """Plot histograms of the psi scores of one label"""
#         fig, axes = plt.subplots(ncols=self.n_clusters, figsize=(4*self.n_clusters,4))
#         for ax, label in zip(axes, self.labels_unique):
        #     print 'label',label
        ax.hist(self.data.psi.ix[self.data.psi.index[self.labels == label],:].values.flat,
                bins=np.arange(0, 1.05, 0.05), facecolor=color, linewidth=0.1)
        ax.set_title('Cluster: {}'.format(label))
        ax.set_xlim(0,1)
        sns.despine()

    def lavalamp(self, ax, label, color):
        """makes a lavalamp of one label"""
        nrow = (self.labels == label).sum()
        ncol = self.data.psi.shape[1]
#         print 'nrow, ncol', nrow, ncol
        x = np.vstack(np.arange(nrow) for _ in range(ncol))
#             x = x+ x_prev
        y = self.data.psi.ix[self.data.psi.index[self.labels == label],:].values.T
#         switchy_scores = np.apply_along_axis(switchy_score, axis=0, arr=y)
#         order = np.argsort(switchy_scores)
        order = get_switchy_score_order(y)
#             x = x[:,order]
        y = y[:,order]
        x_prev = x.max() + 1
        ax.scatter(x, y, color=color, alpha=0.5, edgecolor='#262626', linewidth=0.1)
        sns.despine()
        ax.set_xlim(0, x_prev)
        ax.set_ylim(0, 1)
        ax.set_title('n = {}'.format(nrow))

    def _annotate_centers(self, ax):
        if type(self.clusterer) == KMeans or type(self.clusterer == FuzzyCMeans):
            # Plot the centroids as a white X
            centroids = self.clusterer.cluster_centers_
#         elif type(self.clusterer) == AffinityPropagation:
#             centroids = self.data.binned[self.clusterer.cluster_centers_indices_,:]
#         elif type(self.clusterer) == DBSCAN:
#             return
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
        # Reset the color cycle in case we already cycled through it
        self.color_cycle = cycle(self.colors)

        for label, color in zip(self.labels_unique, self.color_cycle):
            if label % 10 == 0:
                print 'plotting cluster {} of {}'.format(label, self.n_clusters)
            if label == -1:
                color = 'k'
            n_samples_in_cluster = (self.labels == label).sum()
            if n_samples_in_cluster <= 5:
                continue

            fig = plt.figure(figsize=(16, 4))
            hist_ax = plt.subplot2grid((1, 5), (0,0), colspan=1, rowspan=1)
            lavalamp_ax = plt.subplot2grid((1,5), (0, 1), colspan=4, rowspan=4)

            self.hist(hist_ax, label, color=color)
            self.lavalamp(lavalamp_ax, label, color=color)

    def pca_viz(self):
        """Visualizes the clusters on the PCA of the data (from self.data.pca_binned)"""
        # Step size of the mesh. Decrease to increase the quality of the VQ.
#         h = .02     # point in the mesh [x_min, m_max]x[y_min, y_max].

        # Plot the decision boundary. For that, we will assign a color to each
        x_min, x_max = self.reduced[:, 0].min(), self.reduced[:, 0].max()
        y_min, y_max = self.reduced[:, 1].min(), self.reduced[:, 1].max()
#         xx, yy = np.meshgrid (np.arange (x_min, x_max, h), np.arange (y_min, y_max, h))

#         # Obtain labels for each point in mesh. Use last trained model.
#         Z = self.clusterer.predict(np.c_[xx.ravel(), yy.ravel()])

#         # Put the result into a color plot
#         Z = Z.reshape(xx.shape)
#         print Z
    #     pl.figure (1)
    #     pl.clf ()
        fig, ax = plt.subplots(figsize=(12,8))
#         ax.imshow (Z, interpolation='nearest',
#                   extent=(xx.min(), xx.max(), yy.min(), yy.max()),
#                   cmap=mpl.cm.Set1,
#                   aspect='auto', origin='lower', alpha=0.5)

#         ax.plot(self.data.reduced_binned[:, 0], self.data.reduced_binned[:, 1], 'k.',
#                 markersize=2)

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

    def violinplot_random_cluster_members(self, n=20):
        self.color_cycle = cycle(self.colors)
        for label, color in zip(self.labels_unique, self.color_cycle):
            if label == -1:
                color = 'k'
            these_labels = self.labels == label
            events = np.random.choice(self.data.psi.index[these_labels], size=n)
            y = self.data.psi.ix[events,:].values.T
            order = get_switchy_score_order(y)
            events = events[order]
            fig, axes = plt.subplots(ncols=n, figsize=(n, 2), sharey=True)
            for event, ax in zip(events, axes):
        #         if i % 20 == 0:
                sns.violinplot(self.data.psi.ix[event], bw=0.1, inner='points',
                               color=color, linewidth=0, ax=ax, alpha=0.75)
                ax.set_ylim(0,1)
                ax.set_xticks([])
                ax.set_xlabel(label)
                sns.despine()