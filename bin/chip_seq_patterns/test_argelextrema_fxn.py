#!/usr/bin/env python

import sys
import numpy as np
from itertools import tee, izip

def find_local_peak_throughs(x):
    """
    () --> 
    """
    #all entries in the 1d array a smaller than their neighbors
    #from http://stackoverflow.com/questions/4624970/finding-local-maxima-minima-with-numpy-in-a-1d-numpy-array
    #numpy.r_[True, a[1:] < a[:-1]] & numpy.r_[a[:-1] < a[1:], True]
    from scipy.signal import argrelextrema
    #from scipy.signal import argrelextrema    
    # for local maxima
    maxima = argrelextrema(np.array(x), np.greater)
    minima= argrelextrema(np.array(x), np.less)
    #print max
    #print 'max is',type(max)
    #print 'max(0) is', type(max[0])
    #print 'length of max(0) is', len(max[0])
    # for local minima
    if len(maxima[0]) ==0:
        maxima = None
    #print type(min)
    if len(minima[0]) ==0:
        minima = None  
    #if (len(max) ==2 or len(min))  and type(max[0]==)
    #note, these are the indices of x that are local max / min, try    
    return maxima, minima




vector= [1,1,2,3,4,5,3,3,3,1,2,3]
vector2= [1,1,1,1,2,2,2,2,1,1,1]
vector3 =[1,1,1,1,2,22,2,2,1,1,1] 
#print find_local_peak_throughs(vector)
#print find_local_peak_throughs(vector)[0]
#print find_local_peak_throughs(vector)[1]
#print type(find_local_peak_throughs(vector))
print 'vector1 - with max and min'
print vector
print '3x0'+str(find_local_peak_throughs(vector)[0][0][0])
print '3x0'+str(type(find_local_peak_throughs(vector)[0][0][0]))
print '1x0'+str(find_local_peak_throughs(vector)[1])
print '1x0'+str(type(find_local_peak_throughs(vector)[1]))
print 'vector2 - no max no min'
print vector2
print '1x0'+str(find_local_peak_throughs(vector2)[0])
print '1x0'+str(type(find_local_peak_throughs(vector2)[0]))
print '1x1'+str(find_local_peak_throughs(vector2)[1])
print '1x1'+str(type(find_local_peak_throughs(vector2)[1]))
print 'vector3 - max but no min'
print vector3
print '1x0'+str(find_local_peak_throughs(vector3)[0])
print '1x0'+str(type(find_local_peak_throughs(vector3)[0]))
print '2x0'+str(find_local_peak_throughs(vector3)[0][0])
print '2x0'+str(type(find_local_peak_throughs(vector3)[0][0]))


#print '2x0'+str(find_local_peak_throughs(vector2)[0][0])
#print '2x0'+str(type(find_local_peak_throughs(vector2)[0][0]))

#print '3x0'+find_local_peak_throughs(vector2)[0][0][0]
#print '3x0'+type(find_local_peak_throughs(vector2)[0][0][0])


