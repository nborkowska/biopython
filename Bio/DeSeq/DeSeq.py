import numpy as np

from pandas import *
from scipy import stats as st


def dnbinom(x, size, mean):
    """
    probability mass function 
    with alternative parametrization
    """
    prob = size/(size+mean)
    return st.nbinom.pmf(x, size, prob)


class DSet(object):
    """
    Container for a data
    
    Parameters: 
    
    data: pandas DataFrame object, ex:
    
    >>> data = read_table('./your_data.tab', index_col=0)

    conds: array-like, describing experimental conditions, ex:

    >>> conds = ['A','A','B','B']

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

    DISP_MODES = (
            'max',
            'fit-only',
            'gene-only')

    DISP_METHODS = (
            'pooled',
            'per-condition',
            'blind')
    
    def __init__(self, data, conds, sizeFactors=None):
        
        if isinstance(data, DataFrame):
            """ set experimental conditions as hierarchical index """
            index = MultiIndex.from_tuples(zip(conds, data.columns), \
                    names=['conditions','replicates'])
            new = data.reindex(columns = index)
            for idx in index:
                new[idx] = data[idx[1]]
            
            self.data = new
            self.conds = conds
            self.sizeFactors = sizeFactors
        else:
            raise TypeError("Data must be a pandas DataFrame object!")
    
    def setSizeFactors(self, function=np.median):
        """ params: 
            function - use specific function when estimating the factors,
                       median is the default """
        array = self.data.values
        geometricMean = st.gmean(array, axis=1)
        divided = np.divide(np.delete(array, np.where(geometricMean == 0), \
                axis=0).T, [x for x in geometricMean if x != 0])
        self.sizeFactors = function(divided, axis=1)

    def getNormalizedCounts(self):
        return np.divide(self.values, self.sizeFactors)
        
    def setDispersions(self, method='per-condition', mode='max'):
        """ Set dispersion estimates in a data frame """
        
        if mode not in self.DISP_MODES:
            raise ValueError("Invalid mode. Choose from %s, %s, %s." \
                    % self.DISP_MODES)
        if method not in self.DISP_METHODS:
            raise ValueError("Invalid method. Choose from %s, %s, %s." \
                    % self.DISP_METHODS)
        if self.sizeFactors is None:
            raise ValueError("No size factors available. \
                    Call 'setSizeFactors' first.")

        if method == 'pooled':
            rep = self.data.groupby(lambda x: self.conds.count(x), \
                    level=0, axis=1)
            df = len(set(self.conds)) - len(rep.groups.get(1, []))
            normalized = self.getNormalizedCounts()
            bMean = np.mean(normalized)
            print rep.groups[2]
            pass
        elif method == 'per-condition':
            pass
        else:
            pass
