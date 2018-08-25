#!/usr/bin/env python

# remove randomly X % of all reads

import itertools
import HTSeq
import sys
from collections import defaultdict
import gzip
import numpy
from Bio.Seq import Seq

# go through all reads and count them all 
def count_total_reads(fastq_file):
    """
    (file) -> (int)
    goes through fastq file and returns number of reads present
    """
    total = 0
    for line in fastq_file: 
        if line.startswith('@'):
            total+=1 
    return total
    
# calculate number of reads to remove that corresponds to 60%
def calc_no_reads_2_remove(pc_reads, total_no_reads):
    """
    (int), (int) -> (int)
    calculates number of reads that should be removed
    """
    x= int(total_no_reads) * int(pc_reads) /100.0
    return x
    
if __name__ == "__main__":  
    
    #Input files (
    raw_read_filename1 = sys.argv[1]    
    raw_read_filename2 = sys.argv[2]
    pc_reads_2_remove = int(sys.argv[3])

    in1 = iter( HTSeq.FastqReader(raw_read_filename1))
    in2 = iter( HTSeq.FastqReader(raw_read_filename2))

    input1_count_file = open(sys.argv[1], "r")
    
    #Output files
    output_filename1 = str(pc_reads_2_remove)+"_"+"dr_pc_rem_"+sys.argv[1].split('.gz')[0]
    output_filename2 =  str(pc_reads_2_remove)+"_"+"dr_pc_rem_"+sys.argv[2].split('.gz')[0]
    
    out1 = open(output_filename1, 'w')
    out2 = open(output_filename2, 'w')
    
    #Count total number of reads
    total_no_reads = int(count_total_reads(input1_count_file))
    input1_count_file.close()
    print total_no_reads
    
    #Determine number of reads to remove randomly 
    no_reads_2_remove= calc_no_reads_2_remove(pc_reads_2_remove, total_no_reads)
    
    # randomly generate vector corresponding to 60% removed vector 
    reads_to_remove_list1 = numpy.random.randint(0, total_no_reads, no_reads_2_remove)
    reads_to_remove_list2 = numpy.random.randint(0, total_no_reads, no_reads_2_remove)
    # go through each pair of reads and remove the randomly selected reads 
    read_number = 0
    while True:
        for read1, read2 in itertools.izip(in1, in2):
            if len(read1) == 0:
                print 'EOF!'
                break # EOF
            if read_number not in reads_to_remove_list1:
                out1.write('@'+read1.name +  "\n" + read1.seq + "\n"+ "+ \n" + read1.qualstr + "\n")
            if read_number not in reads_to_remove_list2:  
                out2.write('@'+read2.name + "\n" +  read2.seq + "\n"+ "+ \n" + read2.qualstr + "\n")
            read_number +=1
            
