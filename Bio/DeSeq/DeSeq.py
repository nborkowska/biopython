import numpy as np

from pandas import *
from scipy import stats as st
import statsmodels.api as sm


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
            self.disps = None
        else:
            raise TypeError("Data must be a pandas DataFrame object!")
 
    def setSizeFactors(self, function=np.median):
        """ params: 
            function - use specific function when estimating the factors,
                       median is the default """
        array = self.data.values
        geometricMean = st.gmean(array, axis=1)
        divided = np.divide(np.delete(array, np.where(geometricMean == 0),
            axis=0).T, [x for x in geometricMean if x != 0])
        self.sizeFactors = Series(function(divided, axis=1),
                index=self.data.columns)
         
    @staticmethod
    def getNormalizedCounts(dataframe, factors):
        """ factors: array-like or Series """
        return dataframe.div(factors)
    
    @staticmethod
    def getBaseMeansAndVariances(dataframe, factors):
        return DataFrame({
            'bMean': np.mean(dataframe.values, axis=1),
            'bVar': np.var(dataframe.values, axis=1, ddof=1)
            }, index=dataframe.index)

    def selectReplicated(self, normalized):
        return normalized.select(lambda x:self.conds.count(x[0]) > 1,
                axis=1).groupby(axis=1, level=0)

    def _estimateAndFitDispersions(self, mav, sizeFactors, xim):
        dispsAll = Series(
                (mav.bVar - xim * mav.bMean)/(mav.bMean)**2,
                index = mav.index)
        toDel = np.where(mav.bMean.values <= 0)
        dataframe = DataFrame({
            'means': np.log(np.delete(mav.bMean.values, toDel)),
            'variances':np.delete(mav.bVar.values, toDel)
            })
        fit = sm.GLM.from_formula(
                formula='variances ~ means',
                df=dataframe,
                family = sm.families.Gamma(link=sm.families.links.log)).fit()
        return dispsAll, fit
    
    def _calculateDispersions(self, mav, sizeFactors, testing, mode):
        xim = np.mean(1/sizeFactors)
        estDisps, fit = self._estimateAndFitDispersions(
                mav, sizeFactors, xim)
        tframe = DataFrame({'means':np.log(testing)})
        fittedDisp= np.clip(
                (fit.predict(tframe)-xim*testing)/testing**2,
                1e-8, float("Inf"))
        if mode == 'max':
            return Series(np.maximum(estDisps, fittedDisp))
        elif mode == 'fit-only':
            return Series(fittedDisp)
        else:
            return Series(estDisp)

    def setDispersions(self, method='per-condition', mode='fit-only'):
        """ Get dispersion estimates """
        
        if mode not in self.DISP_MODES:
            raise ValueError("Invalid mode. Choose from %s, %s, %s." \
                    % self.DISP_MODES)
        if method not in self.DISP_METHODS:
            raise ValueError("Invalid method. Choose from %s, %s, %s." \
                    % self.DISP_METHODS)
        if self.sizeFactors is None:
            raise ValueError("No size factors available. \
                    Call 'setSizeFactors' first.")

        normalized = DSet.getNormalizedCounts(
                self.data, self.sizeFactors)
        overallBMeans = np.mean(normalized.values, axis=1)
        dfr = {}

        if method == 'pooled':
            """ select all conditions with replicates and estimate a
            single pooled empirical dispersion value """
        
            replicated = self.selectReplicated(normalized)
            groupv = replicated.agg(lambda x: sum((x - np.mean(x))**2))
            bVar = groupv.sum(axis=1) / len(replicated.groups)
            meansAndVars = DataFrame({'bMean':overallBMeans,'bVar':bVar},
                    index=self.data.index) 
            dispersions = self._calculateDispersions(
                    meansAndVars, self.sizeFactors,
                    overallBMeans, mode)
            for name, df in self.data.groupby(axis=1, level=0):
                dfr[name] = dispersions
        
        elif method == 'per-condition':
            replicated = self.selectReplicated(normalized)
            if not replicated.groups:
                raise Exception("None of your conditions is replicated." 
                        + " Use method='blind' to estimate across conditions")
            for name, df in replicated:
                sizeFactors = self.sizeFactors[name].values
                meansAndVars = DSet.getBaseMeansAndVariances(
                        df, sizeFactors)
                dispersions = self._calculateDispersions(
                        meansAndVars, sizeFactors,
                        overallBMeans, mode)
                dfr[name] = dispersions
        else:
            meansAndVars = DSet.getBaseMeansAndVariances(
                    self.data, self.sizeFactors)
            dispersions = self._calculateDispersions(
                    meansAndVars, self.sizeFactors,
                    overallBMeans, mode)
            for name, df in self.data.groupby(axis=1, level=0):
                dfr[name] = dispersions

        self.disps = DataFrame(dfr)
