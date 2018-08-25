#!/usr/bin/env python

import sys
import numpy as np
from itertools import tee, izip

def find_1st_line():
    """(file)-> (str)
    track type=wiggle_0 name=H3K4me3_dome maxHeightPixels=64:64:11 color=31,120,180 visibility=full
    """
    wigfile = open(sys.argv[1],'r')
    for line in wigfile: 
        if line.startswith('track'):
            first_line = line 
            break
    wigfile.close()
    return line
    
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
    #dictionary = {}
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
    minima= argrelextrema(np.array(x), np.less)
    #print max
    #print 'max is',type(max)
    #print 'max(0) is', type(max[0])
    #print 'length of max is', len(max)
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

def parse_wig_decl_track(title):
    """
    (str)--> (str),(int),(int)
    Parses the chromosome number, the position of the first base
    in the peak
    fixedStep chrom=chr1 start=10781 step=10 span=10
    """
    step_type = title.partition(" ")[0]
    chrom = title.split(" ")[1].split('=')[1]
    start = title.partition(" ")[2].partition('=')[1]
    return step_type, chrom, start

def write_new_decl_track(old_decl_track, new_start):
    """
    (string)-> (string)
    fixedStep chrom=chr1 start=10781 step=10 span=10
    """    
    step_type, chrom, old_start = parse_wig_decl_track(old_decl_track)
    new_decl_track = str(str(step_type)+ ' chrom='+str(chrom)+' start='+str(new_start)+' step=10'+' span=10')
    return new_decl_track

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


if __name__ == "__main__":
    from scipy.signal import argrelextrema
    wigfilename= sys.argv[1]    #name of wig input file to parse
    out_wig_filename = 'Better_wig_'+wigfilename.split('.')[0]+'.wig'
    out_bed_filename = 'Better_bed_'+wigfilename.split('.')[0]+'.bed'
    out_out_filename = 'Betterness_'+wigfilename.split('.')[0]+'.out'
    out_file = open(out_out_filename, 'w')
    wig = open(wigfilename,'r') #wig input handler
    out_wig = open(out_wig_filename, 'w')
    out_bed = open(out_bed_filename, 'w')
    out_wig.write(find_1st_line())
    out_bed.write('track name='+wigfilename.split('.')[0]+'\n')
    no_max_no_min = 0
    multi_max = 0
    multi_min = 0 
    no_peaks = 0
    for declaration_track, peak_intensity in read_wig(wig): #for each old 'peak'
        no_peaks += 1 
        chrom = declaration_track.split(" ")[1].split("=")[1]
        len_peak = len(peak_intensity) 
        old_start = declaration_track.split(" ")[2].split("=")[1]
        if (find_local_peak_throughs(peak_intensity)[0] == None) and (find_local_peak_throughs(peak_intensity)[1] == None): #if there are no local maxima or minima in the "peak"
            #x = range(start, start+len_peak, len_peak+1)
            no_max_no_min +=1 
            peak_int = max(peak_intensity)  #use maximum intensity in the whole peak as intensity value
            new_start = ((int(old_start)*2)+(int(len_peak)*10))/2   #new start pos = average position in peak less peak or old_start+end(=old_start+len_peak*10)/2
            new_end = int(new_start)+10 #span and step of 10 in original wig
            out_wig.write(write_new_decl_track(declaration_track, new_start)+'\n'+str(peak_int)+'\n') # write declaration line /n peak intensity /n
            out_bed.write(str(chrom)+'\t'+str(new_start)+'\t'+str(new_end)+'\n')
        elif len(argrelextrema(np.array(peak_intensity), np.greater)[0].tolist()) > 1:       #if there is a maximum or multiple maxima
            multi_max+=1
            maxima = argrelextrema(np.array(peak_intensity), np.greater)[0].tolist()
            for maximum in maxima:
                       new_start = int(declaration_track.split(" ")[2].split("=")[1])+(int(maximum)*10)
                       new_end = int(new_start)+10
                       peak_int= peak_intensity[maximum]
                       out_wig.write(write_new_decl_track(declaration_track, new_start)+'\n'+str(peak_int)+'\n') # write declaration line /n peak intensities /n
                       out_bed.write(str(chrom)+'\t'+str(new_start)+'\t'+str(new_end)+'\n')
        elif type(find_local_peak_throughs(peak_intensity)[1] == tuple) and (find_local_peak_throughs(peak_intensity)[0] == None): #if multiple minima but no maxima
            multi_min +=1
            for v, w in pairwise(find_local_peak_throughs(peak_intensity)[1]):
                new_start = old_start + peak_intensity[v]*10
                new_end  = new_start + peak_intensity[w]*10
                peak_int =  max(peak_intensity)
                out_wig.write(write_new_decl_track(declaration_track, new_start)+'\n'+str(peak_int)+'\n') # write peak
                out_bed.write('Chr'+'\t'+ chrom + '\t' + new_start + '\t' + new_end + '\n')
    
    out_file.write('Went through file: '+ wigfilename + '\n')
    out_file.write('Found maxima in:'+ str(multi_max)+'\n')
    out_file.write('Found minima in:'+ str(multi_min) + '\n')
    out_file.write('Found no maxima and no minima in:'+ str(no_max_no_min) + '\n')
    out_file.write('Total peaks:'+ str(no_peaks) + '\n')