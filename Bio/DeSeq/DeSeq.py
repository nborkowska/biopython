import numpy as np

from pandas import *
from scipy import stats as st
import spam

class DESeq(object):
    """
    Parameters: 
    
    data: pandas DataFrame object, ex:
    
    >>> data = read_table('./your_data.tab', index_col=0)

    conditions: pandas Series object describing experimental conditions, ex:

    >>> conds = Series(["a","a","b","b"], data.columns)

    You may need to add row names to your data:
    
    >>> data.index = ['gene_%d' % x for x in xrange(len(data.values))]
    >>> print data
            untreated1     untreated2     treated1     treated2
    gene_0          56             44           11           23
    gene_1         345            424          560          675
    gene_2          12             45           32           17

    Or create them by mixing some columns:

    >>> print data
    chrom     start     stop    10847_2    10847_3    10847_4
    chr 1    713615   714507         38         75        390
    chr 1    742153   742162         58         11         34
    
    >>> data.index = data.pop('chrom')+':'+data.pop('start').map(str) \ 
    ...     +'-'+data.pop('stop').map(str)
    >>> print data
                           10847_2    10847_3    10847_4
    chr 1:713615-714507         38         75        390
    chr 1:742153-742162         58         11         34 
    """

    def __init__(self, data, conditions, sizeFactors=None):
        
        if isinstance(data, DataFrame):
            self.data = data
            self.conditions = conditions
            self.sizeFactors = sizeFactors
        else:
            raise TypeError("Data must be a pandas DataFrame object!")

    def setSizeFactors(self, function=np.median):
        """ params: 
            function - use specific function when estimating the factors
        """
        array = self.data.values
        geometricMean = st.gmean(array, axis=1)
        divided = np.divide(np.delete(array, np.where(geometricMean == 0), \
                axis=0).T, [x for x in geometricMean if x != 0])
        self.sizeFactors = dict(zip(self.data.columns, \
                function(divided, axis=1)))
