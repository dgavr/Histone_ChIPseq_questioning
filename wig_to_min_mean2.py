#!/usr/bin/env python

import sys
import numpy as np
from itertools import tee, izip

#Mean!

def find_1st_line():
    """(file)-> (str)
    track type=wiggle_0 name=H3K4me3_dome maxHeightPixels=64:64:11 color=31,120,180 visibility=full
    """
    wigfile = open(sys.argv[1],'r')
    for line in wigfile: 
        if line.startswith('track'):
            first_line = 'track type=wiggle_0 name=Minus_mean_'+line.split()[2].split('=')[1]+' maxHeightPixels=64:64:11 color=31,120,180 visibility=full \n'
            break
    wigfile.close()
    return first_line
    
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

def over_threshold(mylist, threshold):
    """
    """
    good = [index for index,value in enumerate(mylist) if value > threshold]
    return good
            
if __name__ == "__main__":
    from scipy.signal import argrelextrema
    from operator import itemgetter
    from itertools import groupby
    
    wigfilename= sys.argv[1]    #name of wig input file to parse
    out_wig_filename = 'Minus_mean_wig_'+wigfilename.split('.')[0]+'.wig'
    out_bed_filename = 'Minus_mean_'+wigfilename.split('.')[0]+'.bed'
    wig = open(wigfilename,'r') #wig input handler
    out_wig = open(out_wig_filename, 'w')
    out_bed = open(out_bed_filename, 'w')
    out_wig.write(find_1st_line())
    out_bed.write('track name='+'Minus_mean' + wigfilename.split('.')[0] + '\n')
    for declaration_track, peak_intensity in read_wig(wig): #for each old 'peak'
        #print 'ORIGINAL TRACK:'+declaration_track +'\n'
        #print 'Original peaks:  \n'
        #print peak_intensity
        threshold = np.mean(peak_intensity)
        chrom = declaration_track.split(" ")[1].split("=")[1]
        old_start = declaration_track.split(" ")[2].split("=")[1]
        start = old_start
        #print 'Mean is: ' +str(threshold)+ '\n'
        new = over_threshold(peak_intensity, threshold)
        #print new
        peak_ranges = []
        for k, g in groupby(enumerate(new), lambda (i,x):i-x):
            group = map(itemgetter(1), g)
            peak_ranges.append((group[0], group[-1]))
        #print peak_ranges
        for peaks in peak_ranges:
            new_start_index= peaks[0]
            new_end_index = peaks[1]
            new_start= int(old_start)+int(new_start_index*10)
            new_end = int(old_start)+int(new_end_index*10)
            #print 'peaks', peaks
            #print 'peaks[0]', peaks[0]
            #print 'peaks[1]', peaks[1]
            new_peak_intensity = peak_intensity[new_start_index:new_end_index]
            out_wig.write(str(write_new_decl_track(declaration_track, new_start))+'\n')
            for item in new_peak_intensity:
              out_wig.write("%s\n" % item)
            out_bed.write(str(chrom)+'\t'+ str(new_start)+'\t'+str(new_end)+'\n')
            
        
        
        