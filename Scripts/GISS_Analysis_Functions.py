import glob, os, re  
from cdo import * 
from netCDF4 import Dataset
from pylab import * 
#import pylab
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.pyplot import * 
from mpl_toolkits.basemap import Basemap



"""
Functions for GISS_Analysis.def foo():
    doc = "The foo property."
    def fget(self):
        return self._foo
    def fset(self, value):
        self._foo = value
    def fdel(self):
        del self._foo
    return locals()
foo = property(**foo())
"""

def common_prefix(strings):
    """ Find the longest string that is a prefix of all the strings.
    """
    if not strings:
        return ''
    prefix = strings[0]
    for s in strings:
        if len(s) < len(prefix):
            prefix = prefix[:len(s)]
        if not prefix:
            return ''
        for i in range(len(prefix)):
            if prefix[i] != s[i]:
                prefix = prefix[:i]
                break
    return prefix






