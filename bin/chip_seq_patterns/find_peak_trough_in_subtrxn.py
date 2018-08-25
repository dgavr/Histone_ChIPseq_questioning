#!/usr/bin/env python
import sys
import numpy 
from collections import defaultdict


def get_chrom_size():
    """
    """
    chrom_dict = defaultdict(list)
    for line in open('chromInfo.txt'):
        chrom = line.split('\t')[0]
        size = line.split('\t')[1]
        chrom_dict[chrom]=size
    return chrom_dict
           
def get_ws_lines_from_file(file, window_size, chrom_dict):
    """
    (file) -> (strings)
    """    
    chrom = file.split('_')[1]
    size_chrom = chrom_dict[chrom]
    nb_it= int(size_chrom)/ int(window_size)
    remainder =  int(size_chrom) % int(window_size)
    line_counter =0
    it_counter =0 
    me1 = []
    ac27 = []
    for line in open(file):
        if line.startswith('Chromosome'):
            print 'skipping title line '
            continue
        else: 
            line_counter+=1
            it_counter += 1
            me1.append(int(line.split('\t')[2]))
            ac27.append(int(line.split('\t')[3]))
            if len(me1) == int(window_size): 
                yield me1, ac27, window_size
                line_counter = 0 
            elif it_counter == nb_it: 
                window_size = remainder-1
            
if __name__ == "__main__":
    subtracted_filename = sys.argv[1]
    print 'opened the file', subtracted_filename
    chrom = subtracted_filename.split('_')[1]
    print 'it is on chromosome', chrom 
    chrom_dict = get_chrom_size()
    print 'finished getting chromosome sizes'
    print chrom, 'contains', str(chrom_dict[chrom]), 'bps'
    window_size = sys.argv[2]
    print 'the window size is ', str(window_size)
    output_filename = 'peaks_'+subtracted_filename[subtracted_filename.find('sub')+4:]+'.bed'
    output_file = open(output_filename, 'w')
    print 'finished creating output file called', output_filename 
    output_file.write('track name='+output_filename.split('.')[0]+'\n')
    window_counter= 1
    start_pos = 0
    for me1, ac27, new_window_size in get_ws_lines_from_file(subtracted_filename, window_size, chrom_dict):
        print 'yielded first', str(window_size), 'bps'
        print 'Me1 list from def is:', me1
        print 'Ac27 list fron def is:', ac27
        me1_indices =  numpy.diff(numpy.sign(numpy.diff(me1)))
        print 'finished calculating indices of me1'
        print 'me1 indices looks like this', me1_indices
        ac27_indices = numpy.diff(numpy.sign(numpy.diff(ac27)))
        print 'finished calculating indices of ac27'
        print 'ac27 indices looks like this', ac27_indices
        peaks_me1 = numpy.where(me1_indices < -1)
        print 'Me1 Peaks is a', type(peaks_me1), 'and looks like this', peaks_me1
        peaks_me1 = peaks_me1[0]
        peaks_ac27 = numpy.where(ac27_indices< -1) 
        print 'Ac27 Peaks is a', type(peaks_me1), 'and looks like this', peaks_ac27
        peaks_ac27= peaks_ac27[0]
        print 'Ac 27[0] looks like this', peaks_ac27
        print 'Me1 peaks[0] is a ', type(peaks_me1), 'and looks like this', peaks_me1
        #overlap= overlap(peaks_me1, peaks_ac27)
        common_peaks = numpy.intersect1d(peaks_me1, peaks_ac27)
        print 'Common peaks are', type(common_peaks), 'and look like this:', common_peaks
        print 'Set of common_peaks: ', set(common_peaks)
        print 'sorted list of set of common_peaks', sorted(list(set(common_peaks)))
        #common_peaks = set(peaks_me1.flatten).intersection(set(peaks_ac27.flatten))
        for common_peak in common_peaks: 
                pos = int(common_peak) + int(window_counter)*int(window_size)
                output_file.write('\t'.join([str(chrom), str(pos+1), str(pos+2)])+'\n')
        window_counter+=1    
