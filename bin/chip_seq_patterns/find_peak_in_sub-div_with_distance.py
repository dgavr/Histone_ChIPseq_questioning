#!/usr/bin/env python
import sys
import numpy 
from collections import defaultdict


def zero_if_lower_than0(x):
    if x < 0:
        return 0
    return x


def get_chr_size(chrom_size_filename):
    size_dictionary = defaultdict(int)
    for line in open(chrom_size_filename):
        chrom ,size_of_chrom =line.rstrip().split('\t')
        #print chrom, size_of_chrom
        size_dictionary[chrom] = int(size_of_chrom)
    return size_dictionary

           

def go_thru_bedgraph(bedgraph_filename):
    bg_file = open(bedgraph_filename,'r')
    bg_dictionary = defaultdict(list)
    for i,line in enumerate(bg_file):
        if not line.startswith('track'):
            chrom,start_pos,end_pos,intensity = line.strip().split('\t')
            bg_dictionary[chrom].append((int(float(start_pos)),int(float(end_pos)), int(float(intensity))))
            #print bg_dictionary[chrom]
    return bg_dictionary
    
def compare_chr(dictionary1, dictionary2):
    """ Make chromosomes (the keys to the dictionary generated in go_thru() and export
        them as a list. Then compare the list of chromosomes to make sure that they are 
        all present in both bed files. Export the list of K present in both files. Gets 
        rid of chromosomes absent in one of the files.
    """
    intersect = set(dictionary1.keys()).intersection(set(dictionary2.keys()))
    chr_list = list(intersect)
    return chr_list


def fill(size_of_chrom,bedlist):
    list_pos_in_chrom = [0]*size_of_chrom
    for peak in bedlist:
        for pos in range(peak[0],peak[1]):
            list_pos_in_chrom[pos] = peak[2]
    return list_pos_in_chrom


def find_peak_within_dist(list1, list2, distance):
    """
    """
    new_list=[]
    #no_pairs = 0
    i = 0
    a = 0
    b = 0
    if distance < len(list1):
        while i < len(list1):
            print list1[i], type(list1[i])
            print list2[a], type(list2[a])
            print list2[b], type(list2[b])
            while ((a < i) and (list1[i]-list2[a] > int(distance))):
                #new_list.append((str(list1[i]), str(list2[a])))
                a +=1
            while ((b < len(list2)) and (list2[b]-list1[i]< int(distance))):
                print 'second criteria true'
                new_list.append((str(list2[b]), str(list1[i])))
                b+=1
            i +=1
    return new_list
    
                
if __name__ == "__main__":
    #Input filenames
    histone_mark1_filename = sys.argv[1]
    histone_mark2_filename = sys.argv[2]
    chrom_size_filename = sys.argv[3]
    distance= int(sys.argv[4])
    print 'opened the files:', histone_mark1_filename, histone_mark2_filename
    
    print "I'm going to make dictionaries with them now "
    print 'Now loading into dictionary the file: ', histone_mark1_filename
    
    bed1_dictionary = go_thru_bedgraph(histone_mark1_filename)
    print 'Finished doing that!'
    print 'Now loading the file ',histone_mark2_filename, 'to dictionary' 
    bed2_dictionary = go_thru_bedgraph(histone_mark2_filename)
    print 'Yay! finished both dictionaries'
    print 'Now making dictionary with chromosome sizes '
    size_dictionary = get_chr_size(chrom_size_filename)
    print 'Finished loading chromosome sizes '
    
    #Make sure the same set of chromosomes is present in both datasets 
    list_of_common_chrom = compare_chr(bed1_dictionary, bed2_dictionary)
    print 'Finished looking for common chromosomes'
    print 'Now going to look for ndrs and peaks!'

    #Make outputfile
    output_peak_filename = 'peak_'+histone_mark1_filename.split('.')[0]+histone_mark2_filename.split('.')[0]+'.txt'
    output_peak_file = open(output_peak_filename, 'w')
    output_peak_file.write('track name='+output_peak_filename.split('.')[0]+ '\n')
    
    window_counter= 1
    start_pos = 0
    for chrom in size_dictionary:
        print 'chrom is', chrom
        #Make lists with the intensities on each chromosome
        histone_mark1_list = fill(size_dictionary[chrom],bed1_dictionary[chrom])
        histone_mark2_list = fill(size_dictionary[chrom],bed2_dictionary[chrom])
        clipped_hm1_list =  [zero_if_lower_than0(x) for x in histone_mark1_list]
        clipped_hm2_list = [zero_if_lower_than0(x) for x in histone_mark2_list]
        hm1_indices =  numpy.diff(numpy.sign(numpy.diff(clipped_hm1_list)))
        print 'finished calculating indices of histone mark1'
        print 'me1 indices looks like this', hm1_indices
        hm2_indices = numpy.diff(numpy.sign(numpy.diff(clipped_hm2_list)))
        print 'finished calculating indices of histone mark2'
        print 'ac27 indices looks like this', hm2_indices
        peaks_hm1 = numpy.where(hm1_indices < -1)
        print 'Histone mark1 peaks is a', type(peaks_hm1), 'and looks like this', peaks_hm1
        peaks_hm1 = peaks_hm1[0]
        peaks_hm2 = numpy.where(hm2_indices< -1) 
        print 'Histone mark2 peaks is a', type(peaks_hm2), 'and looks like this', peaks_hm2
        peaks_hm2= peaks_hm2[0]
        print 'HM2[0] looks like this', peaks_hm2
        print 'HM1 peaks[0] is a ', type(peaks_hm1), 'and looks like this', peaks_hm1
        new_peaks=find_peak_within_dist(peaks_hm1, peaks_hm2, distance)                    
        #overlap= overlap(peaks_me1, peaks_ac27)
        #common_peaks = numpy.intersect1d(peaks_hm1, peaks_hm2)
        #print 'Common peaks are', type(common_peaks), 'and look like this:', common_peaks
        #print 'Set of common_peaks: ', set(common_peaks)
        #print 'sorted list of set of common_peaks', sorted(list(set(common_peaks)))
        #common_peaks = set(peaks_me1.flatten).intersection(set(peaks_ac27.flatten))
        for new_peak in new_peaks: 
                print 'pos is ', type(new_peak)
                output_peak_file.write('\t'.join([str(chrom), str(new_peak[0]), str(new_peak[1])])+'\n')
        window_counter+=1    
