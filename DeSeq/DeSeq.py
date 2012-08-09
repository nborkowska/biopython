import numpy as np

from pandas import *
from scipy import stats as st


class DESeq(object):
    """
    Parameters: 
    
    data: pandas DataFrame object, ex:
    
    >>> from DeSeq import *
    >>> data = read_table('./your_data.tab')
    >>> DESeq(data)

    First column is responsible for row names, ex:
    >>> print data
    gene_id  cond_A  cond_B
      X4222     546      45
      X3445     11       23 
    
    If your data doesn't have such a column, you should add one,
    ex like this:

    >>> data.insert(0,'gene',['gene_%d' % x for x in xrange(len(data.values))])
    
    You may also need to rename some columns:
    
    >>> data = data.rename(columns={'Unnamed: 0':'Genes'})
    
    Or edit some columns like this:

    >>> print data
    chrom     start     stop    10847_2    10847_3    10847_4
    chr 1    713615   714507         38         75        390
    chr 1    742153   742162         58         11         34
    
    >>> new = data.chrom + ":" + data.pop('start').map(str) \
    ...     + "-" + data.pop('stop').map(str)
    >>> data['chrom'] = new
    >>> print data
                  chrom    10847_2    10847_3    10847_4
    chr 1:713615-714507         38         75        390
    chr 1:742153-742162         58         11         34 
    """

    def __init__(self, data, sizeFactors=None):
        
        if isinstance(data, DataFrame):
            self.data = data
            self.sizeFactors = sizeFactors
        else:
            raise TypeError("Data must be a pandas DataFrame object!")

    def setSizeFactors(self, omit=1, function=np.median):
        """ params: 
            omit: omit n first columns corresponding to row names
            function: use specific function when estimating the factors
        """
        array = self.data.values[:,omit:].astype(float)
        geometricMean = st.gmean(array, axis=1)
        divided = np.divide(np.delete(array, np.where(geometricMean == 0), \
                axis=0).T, [x for x in geometricMean if x != 0])
        self.sizeFactors = function(divided, axis=1)
