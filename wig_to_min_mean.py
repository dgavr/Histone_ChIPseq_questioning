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
    from operator import itemgetter
    from itertools import groupby
    
    wigfilename= sys.argv[1]    #name of wig input file to parse
    out_wig_filename = 'Better3_wig_'+wigfilename.split('.')[0]+'.wig'
    out_bed_filename = 'Better3_bed_'+wigfilename.split('.')[0]+'.bed'
    out_out_filename = 'Betterness3_'+wigfilename.split('.')[0]+'.out'
    out_file = open(out_out_filename, 'w')
    wig = open(wigfilename,'r') #wig input handler
    out_wig = open(out_wig_filename, 'w')
    out_bed = open(out_bed_filename, 'w')
    out_wig.write(find_1st_line())
    out_bed.write('track name='+'Try2' + wigfilename.split('.')[0] + '\n')
    no_max_no_min = 0
    multi_max = 0
    multi_min = 0 
    dunno = 0
    for declaration_track, peak_intensity in read_wig(wig): #for each old 'peak'
        print 'ORIGINAL TRACK:'+declaration_track +'\n'
        print 'Original peaks:  \n'+ str(peak_intensity)
        chrom = declaration_track.split(" ")[1].split("=")[1]
        len_peak = len(peak_intensity) 
        old_start = declaration_track.split(" ")[2].split("=")[1]
        peak_ranges = []
        print 'Mean is: ' +str(np.mean(peak_intensity))+ '\n'
        for k, g in groupby(enumerate([x[0] for x in enumerate(peak_intensity) if x[1] <= np.mean(peak_intensity)]), lambda (i,x):i-x):
            print '\n'
            group= map(itemgetter(1),g)
            #print group
            #print 'Group \n'
            peak_ranges.append((group[0], group[1], group[-1]))
            print 'Ranges: '+ str(peak_ranges)
            size= len(peak_ranges)
            print 'Size= '+ str(size)
            beg= 0
            start =old_start
            for p in range (0, size):
                print 'NEW_TRACK:' + str(write_new_decl_track(declaration_track, start)) + '\n'
                peak = peak_intensity[beg:peak_ranges[p][0]]
                print str(peak)+'\n'
                beg = peak_ranges[p][1]
                start = int(start) + int(beg)*10 
                print '\n\n'
                if p == size and p> 1:
                   peak= peak_intensity[beg:-1] 
                   print 'LAST_TRACK:' + str(write_new_decl_track(declaration_track, start))+'\n'
                   print peak
                   print '\n'