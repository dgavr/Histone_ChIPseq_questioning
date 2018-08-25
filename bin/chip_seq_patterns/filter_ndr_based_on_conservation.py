#!/usr/bin/env python
import sys
from collections import defaultdict

def go_thru_bed(bed_filename):
    """(bed file) -> (dictionary)
        Bed obtained using BamToBed from bedtools
        chr1    701     747     WTCHG_45716_234:6:1104:14793:11212#TCCAGCCA/1   3       +
        chr1    703     749     WTCHG_45716_234:6:2201:19323:2626#TCCAGCCA/1    3       +
        chr1    972     1018    WTCHG_45716_234:6:1205:20672:166536#TCCAGCCA/1  40      +
        chr1    973     1019    WTCHG_45716_234:6:1305:17188:58002#TCCAGCCA/1   42      -
        chr1    975     1021    WTCHG_45716_234:6:2206:9505:64053#TCCAGCCA/1    8       +
        dictionary[chr1] = 701-747
    """
    dictionary = defaultdict(list) 
    bedfile = open(bed_filename,'r')
    for i,line in enumerate(bedfile):
        if not line.startswith('track'):
           chrom = line.rstrip().split('\t')[0]
           start_pos =line.rstrip().split('\t')[1]
           end_pos = line.rstrip().split('\t')[2]
           dictionary[chrom].append((int(start_pos),int(end_pos)))       
    return dictionary 


def compare_chr(dictionary1, dictionary2):
    """ Make chromosomes (the keys to the dictionary generated in go_thru() and export
        them as a list. Then compare the list of chromosomes to make sure that they are 
        all present in both bed files. Export the list of K present in both files. Gets 
        rid of chromosomes absent in one of the files.
    """
    intersect = set(dictionary1.keys()).intersection(set(dictionary2.keys()))
    chr_list = list(intersect)
    return chr_list


if __name__ == "__main__":
    #Declare inputs
    input_filename = sys.argv[1]
    conservation_filename = sys.argv[2]
        
    #Import input files into dictionaries
    input_dict = go_thru_bed(input_filename)
    cons_dict = go_thru_bed(conservation_filename)
    
    #Determine list of chromosomes common to both datasets
    common_chr = compare_chr(input_dict, cons_dict)
    
    