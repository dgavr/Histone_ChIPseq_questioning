#!/usr/bin/env python

import sys
import numpy as np

def read_wig(wig):
    """
    (file)--> (tuple)
    Opens wig file, reads line and parses
    fixedStep chrom=chr1 start=10781 step=10 span=10 
    """    
    last= None               # buffer for keeping last unprocessed line 
    while True:              # mimic closure??
        if not last:         # the first record or a record following a declaration line
            for line in wig:    # search for start of next peak
#                 if 'track' in line: 
#                    yield parse_wig_track_line(line)
                 if 'fixedStep' in line: # beginning of  declaration track
                    last = line[:-1]     #save this line
                    break   # leaves for loop
        if not last: break  # leaves while loop
        declaration_track, peak_intensity, last = last, [], None   #initializes declaration track name, peak_intensity list and clears last
        for line in wig: # read the sequence
            if  'fixedStep' in line:    #declaration line
                yield declaration_track, peak_intensity
                last = line[:-1]            #save this line
                break               #leave for loop
            peak_intensity.append(float(line[:-1]))#no else! if it's not a declaration line then it's a peak intensity value, stays in for loop until reaches a line with 'fixedStep' which is a declaration track

def count_unsorted_list_items(items):
    """
    (tuple)--> ()
    e.g. tuple = (1,2,1,1,1,10,2,3,4,10)
    Uses tuple from read_wig() fuction and counts the number of occurences of each peak intensity within peak
    Code from http://stackoverflow.com/questions/2600191/how-to-count-the-occurrences-of-a-list-item-in-python
    """
    #Distribution
    from collections import defaultdict
    dictionary = {}
    counts = defaultdict(int)
    for item in items:
        counts[item] += 1
    return dict(counts)

# def plot_distribution(dictionary_distrib):
#     """
#     (dictionary)-->(plot)
#     Plots the distribution of counts of each peak intensity as a histogram
#     From http://stackoverflow.com/questions/5926061/plot-histogram-in-python
#     """
#     import numpy as np
#     import matplotlib.pyplot as plt
#     
#     
#     pos = np.arange(len(dictionary_distrb.keys()))
#     width = 1.0     # gives histogram aspect to the bar diagram
#     
#     ax = plt.axes()
#     ax.set_xticks(pos + (width / 2))
#     ax.set_xticklabels(dictionary_distrb.keys())
#     
#     plt.bar(pos, dictionary_distrb.values(), width, color='r')
#     plt.show()
    
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
    #print max
    #print 'max is',type(max)
    #print 'max(0) is', type(max[0])
    #print 'length of max is', len(max)
    #print 'length of max(0) is', len(max[0])
    # for local minima
    if len(maxima[0]) ==0:
        maxima = None
    minima= argrelextrema(np.array(x), np.less) 
    #print type(min)
    if len(minima[0]) ==0:
        minima = None
    #if (len(max) ==2 or len(min))  and type(max[0]==)
    #note, these are the indices of x that are local max / min, try    
    #    x[argrelextrema(x, np.greater)[0]]
    return maxima, minima

def parse_wig_decl_track(title):
    """
    (str)--> (str),(str),(int),(int)
    Parses the chromosome number, the position of the first base
    in the peak
    fixedStep chrom=chr1 start=10781 step=10 span=10
    """
    step_type = title.partition(" ")[0]
    chrom = title.partition(" ")[1].partition('=')[1]
    start = title.partition(" ")[2].partition('=')[1]
    step_no = title.partition(" ")[3].partition('=')[1]
    span = title.partition(" ")[4].partition('=')[1]
    return step_type, chrom, start, step_no, span

        
if __name__ == "__main__":

    wigfilename= sys.argv[1]    #name of wig input file to parse
    #new_wigfilename = 'new'+wigfilename #name of output wig file
    out_filename = 'dist'+wigfilename.split('.')[0]+'.txt'
    wig = open(wigfilename,'r') #wig input handler
    out_file = open(out_filename, 'w')
    #new_wig = 
    #declaration_track, peak_intensity = 0, 0    #initialize declaration_track and peak_intensity
    stop = 0
    for declaration_track, peak_intensity in read_wig(wig):
        chrom = declaration_track.split(" ")[3].split("=")[1]
        num = len(peak_intensity)
        start = declaration_track.split(" ")[2].split("=")[1]
        dist = int(stop) - int(start)
        stop = int(start) + num
        out_file.write(chrom+'\t'+str(dist)+'\n')

        #print declaration_track
        #print len(peak_intensity)
        #print count_unsorted_list_items(peak_intensity)
        #print find_local_peak_throughs(peak_intensity)
       