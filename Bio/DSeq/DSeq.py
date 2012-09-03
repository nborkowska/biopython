import numpy
from pandas import *
import scipy.stats
import statsmodels.api as sm
import pp

ppservers = ()
job_server = pp.Server(ppservers=ppservers)

def dnbinom(x, size, mean):
    """
    probability mass function 
    with alternative parametrization
    """
    prob = size/(size+mean)
    return numpy.exp(scipy.stats.nbinom.logpmf(x, size, prob))

def _calculatepValue(inputData):
    index, row, kss, sumDisps, sfSum, mus = inputData
    if all(v == 0 for v in row.values):
        pval=numpy.nan
    else:
        ks = range(int(row.sum())+1)
        """ probability of all possible counts sums with 
        the same total count """
                
        ps = dnbinom(
                ks, 1/sumDisps[0],
                mus*sfSum[0]
                )*dnbinom(row.sum()-ks, 1/sumDisps[1],
                        mus*sfSum[1])

        """ probability of observed count sums """
        pobs = dnbinom(
                kss[0], 1/sumDisps[0],
                mus*sfSum[0]
                )*dnbinom(kss[1], 1/sumDisps[1],
                        mus*sfSum[1])
                
        if kss[0]*sfSum[1] < kss[1]*sfSum[0]:
            number = ps[:int(kss[0]+1)]
        else:
            number = ps[int(kss[0]):]
        pval = numpy.nanmin([1, 2*numpy.nansum(number)/numpy.nansum(ps)]) 
    return pval


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
 
    def setSizeFactors(self, function=numpy.median):
        """ params: 
            function - use specific function when estimating the factors,
                       median is the default """
        array = self.data.values
        geometricMean = scipy.stats.gmean(array, axis=1)
        divided = numpy.divide(numpy.delete(array, numpy.where(geometricMean == 0),
            axis=0).T, [x for x in geometricMean if x != 0])
        self.sizeFactors = Series(function(divided, axis=1),
                index=self.data.columns)
         
    @staticmethod
    def getNormalizedCounts(dataframe, factors):
        """ factors: array-like or Series """
        return dataframe.div(factors)
    
    @staticmethod
    def getBaseMeansAndVariances(dataframe):
        """ dataframe - DataFrame with normalized data """
        return DataFrame({
            'bMean': numpy.mean(dataframe.values, axis=1),
            'bVar': numpy.var(dataframe.values, axis=1, ddof=1)
            }, index=dataframe.index)

    def selectReplicated(self, normalized):
        return normalized.select(lambda x:self.conds.count(x[0]) > 1,
                axis=1).groupby(axis=1, level=0)

    def _estimateAndFitDispersions(self, mav, sizeFactors, xim):
        dispsAll = Series(
                (mav.bVar - xim * mav.bMean)/(mav.bMean)**2,
                index = mav.index)
        toDel = numpy.where(mav.bMean.values <= 0)
        dataframe = DataFrame({
            'means': numpy.log(numpy.delete(mav.bMean.values, toDel)),
            'variances':numpy.delete(mav.bVar.values, toDel)
            })
        fit = sm.GLM.from_formula(
                formula='variances ~ means',
                df=dataframe,
                family = sm.families.Gamma(link=sm.families.links.log)).fit()
        return dispsAll, fit
    
    def _calculateDispersions(self, mav, sizeFactors, testing, mode):
        xim = numpy.mean(1/sizeFactors)
        estDisp, fit = self._estimateAndFitDispersions(
                mav, sizeFactors, xim)
        tframe = DataFrame({'means':numpy.log(testing)})
        fittedDisp= numpy.clip(
                (fit.predict(tframe)-xim*testing)/testing**2,
                1e-8, float("Inf"))
        if mode == 'max':
            disp = numpy.maximum(estDisp, fittedDisp)
        elif mode == 'fit-only':
            disp = fittedDisp
        else:
            disp = estDisp
        return Series(numpy.maximum(disp, 1e-8), index=mav.index)

    def setDispersions(self, method='per-condition', mode='gene-only'):
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
        overallBMeans = numpy.mean(normalized.values, axis=1)
        dfr = {}

        if method == 'pooled':
            """ select all conditions with replicates and estimate a
            single pooled empirical dispersion value """
        
            replicated = self.selectReplicated(normalized)
            groupv = replicated.agg(lambda x: sum((x - numpy.mean(x))**2))
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
                        df)
                dispersions = self._calculateDispersions(
                        meansAndVars, sizeFactors,
                        overallBMeans, mode)
                dfr[name] = dispersions
            maxDisps = DataFrame(dfr).max(axis=1)
            for name, df in normalized:
                if not dfr.has_key(name):
                    dfr[name] = maxDisps
        else:
            meansAndVars = DSet.getBaseMeansAndVariances(
                    self.data)
            dispersions = self._calculateDispersions(
                    meansAndVars, self.sizeFactors,
                    overallBMeans, mode)
            for name, df in self.data.groupby(axis=1, level=0):
                dfr[name] = dispersions

        dfr = DataFrame(dfr)
        self.disps = dfr.fillna(1e-8)
    
    @staticmethod
    def _calculatepValue(inputData):
        print 'ee'
        if all(v == 0 for v in row.values):
            pval=numpy.nan
        else:
            ks = range(int(row.sum())+1)
            """ probability of all possible counts sums with 
            the same total count """
                
            ps = dnbinom(
                    ks, 1/sumDisps.ix[index, 0],
                    mus[index]*sfSum[0]
                    )*dnbinom(row.sum()-ks, 1/sumDisps.ix[index,1],
                            mus[index]*sfSum[1])

            """ probability of observed count sums """
            pobs = dnbinom(
                    kss.ix[index, 0], 1/sumDisps.ix[index, 0],
                    mus[index]*sfSum[0]
                    )*dnbinom(kss.ix[index, 1], 1/sumDisps.ix[index,1],
                            mus[index]*sfSum[1])
                
            if kss.ix[index,0]*sfSum[1] < kss.ix[index, 1]*sfSum[0]:
                number = ps[:int(kss.ix[index,0]+1)]
            else:
                number = ps[int(kss.ix[index,0]):]
            pval = numpy.nanmin([1, 2*numpy.nansum(number)/numpy.nansum(ps)]) 
        return pval

    def _getpValues(self, counts, sizeFactors, disps):
        
        kss = counts.sum(axis=1, level=0).dropna(axis=1,how='all')
        mus = DSet.getNormalizedCounts(counts, sizeFactors).mean(axis=1)
        sumDisps, pvals = {}, []
        for name, col in counts.groupby(level=0, axis=1):
            n = mus*sizeFactors[name].sum()
            fullVars = numpy.maximum(
                    n + disps[name]*numpy.power(mus,2)*numpy.sum(
                        numpy.power(sizeFactors[name].values, 2)),
                    n*(1+1e-8)
                    )
            sumDisps[name] = (fullVars - n) / numpy.power(n, 2)
       
        sumDisps = DataFrame(sumDisps)
        sfSum = sizeFactors.sum(level=0).dropna()
        inputs = []
        for index, row in kss.iterrows():
            inputs.append([index, row, kss.ix[index], sumDisps.ix[index], sfSum, mus[index]])
        inputs=tuple(inputs)
        print job_server.get_ncpus()
        jobs = [job_server.submit(_calculatepValue, (input, ), (dnbinom, ),
        ("numpy", "scipy.stats",  )) for input in inputs]
        
        pvals = [job() for job in jobs]

        return Series(pvals,index=counts.index)
    
    @staticmethod
    def _BenjaminiHochberg(pvals):
        pvals =  pvals.order(na_last=True, ascending=False).dropna()
        l = len(pvals)
        try:
            previous = pvals[0]
            for i,j in enumerate(pvals.iteritems()):
                corrected = min(j[1]*l/(l-i),previous)
                previous = pvals[i] = corrected
        except:
            pass
        return pvals

    def nbinomTest(self, condA, condB):
        if self.disps is None:
            raise ValueError("No dispersion values available."
                    + " Call 'setDispersions' first.")
        if any(cond not in set(self.conds) for cond in [condA, condB]):
            raise ValueError("No such conditions!")
        
        func = lambda x: x[0] in [condA, condB]
        testingConds = self.data.select(func, axis=1)
        sizeFactors = self.sizeFactors.select(func)
        normalizedConds = DSet.getNormalizedCounts(testingConds, sizeFactors)
        meansAndVars = DSet.getBaseMeansAndVariances(normalizedConds)
        dispersions = self.disps.select(lambda x: x in [condA, condB], axis=1)
        p_vals = self._getpValues(testingConds, sizeFactors, dispersions) 
        adjustedPVals = DSet._BenjaminiHochberg(p_vals)
        bmvA = DSet.getBaseMeansAndVariances(normalizedConds[condA])
        bmvB = DSet.getBaseMeansAndVariances(normalizedConds[condB])
        return DataFrame({
            'baseMean': meansAndVars.bMean,
            'baseVar': meansAndVars.bVar,
            'baseMeanA': bmvA.bMean,
            'baseVarA': bmvA.bVar,
            'baseMeanB': bmvB.bMean,
            'baseVarB': bmvB.bVar,
            'pval': p_vals,
            'pvalAdj': adjustedPVals,
            'foldChange': bmvB.bMean / bmvA.bMean,
            'log2FoldChange': numpy.log2( bmvB.bMean / bmvA.bMean)
            }, index=self.data.index)
