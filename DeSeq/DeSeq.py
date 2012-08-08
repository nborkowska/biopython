from pandas import *

class DESeq(object):
    """
    Parameters: 
    
    data: pandas DataFrame object, ex:
    
    >>> from DeSeq import *
    >>> data = read_table('./your_data.tab')
    >>> DESeq(data)

    You may also need to rename some columns:
    
    >>> data = data.rename(columns={'Unnamed: 0':'Genes'})
    
    """
    def __init__(self, data):
        self.data = data
