from singlesail import parsers
from singlesail.data_model import ExpressionData, SplicingData

class Study(object):
    """

    Attributes
    ----------


    Methods
    -------

    """

    def __init__(self, sample_info_filename, expression_matrix_filename,
                 splicing_info_filename):
        """Constructor for Study object containing gene expression and
        alternative splicing data.

        Parameters
        ----------
        sample_info_filename, expression_matrix_filename, splicing_info_filename

        Returns
        -------


        Raises
        ------

        """
        self._sample_info_filename = sample_info_filename
        self._expression_matrix_filename = expression_matrix_filename
        self._splicing_info_filename = splicing_info_filename

        self.sample_info = parsers.read_sample_info(sample_info_filename)
        self.expression = ExpressionData(expression_matrix_filename)
        self.splicing = SplicingData(splicing_info_filename)


    def detect_outliers(self):
        """Detects outlier cells from expression, mapping, and splicing data and labels the outliers as such for future analysis.

        Parameters
        ----------
        self

        Returns
        -------


        Raises
        ------

        """
        raise NotImplementedError

    def jsd(self):
        """Performs Jensen-Shannon Divergence on both splicing and expression data

        Jensen-Shannon divergence is a method of quantifying the amount of
        change in distribution of one measurement (e.g. a splicing event or a
        gene expression) from one celltype to another.
        """
        self.expression.jsd()
        self.splicing.jsd()

    def pca(self):
        """Performs PCA on both expression and splicing data
        """
        self.expression.pca()
        self.splicing.pca()