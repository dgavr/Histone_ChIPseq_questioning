#!/usr/bin/env python

# remove randomly X % of all reads

import itertools
import HTSeq
import sys
from collections import defaultdict
import subprocess
import gzip
import numpy
from Bio.Seq import Seq
import random

# go through all reads and count them all 
def count_total_reads(input_filename1, input_filename2):
    """
    (file) -> (int)
    goes through fastq file and returns number of reads present
    """
    in1 = iter(HTSeq.FastqReader(input_filename1))
    in2 = iter(HTSeq.FastqReader(input_filename2))
    
    #Determine number of reads to remove randomly 
    read_number = 0
    #while True:
    for read1, read2 in itertools.izip(in1, in2):
        #print read_number
        #print read1.name
        read_number +=1
        if len(read1) == 0:
            print 'EOF!'
            break # EOF
    #print 'Finished!'
    print 'total number of reads =', read_number 
    
    in1.close()
    in2.close()
    
    return read_number
        
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

    #Count total number of reads
    #total_no_reads = int(count_total_reads(raw_read_filename1, raw_read_filename2))
    
    #print 'The total number of reads in the original file is ',total_no_reads
    
    line_count=subprocess.Popen(['wc','-l', sys.argv[1]],stdout=subprocess.PIPE)
    (std_out, std_err)= line_count.communicate()
    total_no_reads= int(float(std_out.split(' ')[0])/4.0 )
    print 'The total number of reads in the original file is ',total_no_reads
    #Determine number of reads to remove randomly 
    no_reads_2_remove= calc_no_reads_2_remove(pc_reads_2_remove, total_no_reads)
    print 'I wanna remove ', pc_reads_2_remove, ' \%  of the reads which corresponds to ', no_reads_2_remove
    # randomly generate vector corresponding to 60% removed vector  
    a = numpy.arange(int(total_no_reads-1))
    print 'Vector in range total no reads before shuffling', a
    numpy.random.shuffle(a)
    print 'I just shuffled that same vector', a
    reads_to_remove_list = a[0:(no_reads_2_remove)]
    print 'List of reads to remove:',  reads_to_remove_list, 'which is ', len(reads_to_remove_list), 'reads long'
    reads_to_remove_list.sort()
    print 'sorted version', reads_to_remove_list
    #Open input files
    in1 = iter( HTSeq.FastqReader(raw_read_filename1))
    in2 = iter( HTSeq.FastqReader(raw_read_filename2))
    #input1_count_file = open(sys.argv[1], "r")
    
    #Output files
    #output_filename1 =  str(pc_reads_2_remove) + "_" + "pc_rem_" + sys.argv[1].split('.gz')[0]
    #output_filename2 =  str(pc_reads_2_remove) + "_" + "pc_rem_" + sys.argv[2].split('.gz')[0]
    
    output_filename1 =  str(pc_reads_2_remove)+"_"+"pc_rem_"+sys.argv[1]
    output_filename2 =  str(pc_reads_2_remove)+"_"+"pc_rem_"+sys.argv[2]
    
    out1 = open(output_filename1, 'w')
    out2 = open(output_filename2, 'w')
    
    
    #reads_to_remove_list = numpy.random.randint(0, total_no_reads, no_reads_2_remove)
        # go through each pair of reads and remove the randomly selected reads 
    
    read_number = -1
    vector_counter = 0
    
    #while True:
    for read1, read2 in itertools.izip(in1, in2):
        read_number +=1
        #print 'read_no', read_number
        if (read_number == int(reads_to_remove_list[vector_counter])):
            vector_counter += 1 
            #print 'found one!'
        elif (read_number != int(reads_to_remove_list[vector_counter])):
            out1.write('@'+read1.name +  "\n" + read1.seq + "\n"+ "+ \n" + read1.qualstr + "\n")
            out2.write('@'+read2.name + "\n" +  read2.seq + "\n"+ "+ \n" + read2.qualstr + "\n")
            #print 'I will ignore this read it is read no', read_number, ' and corresponds to ', read1.name
        if len(read1) == 0:
            print 'EOF!'
            break # EOF
    #print 'Finished!'
    print 'Number of reads in original file = ', read_number
    print 'After removing ', pc_reads_2_remove, '% of reads, which corresponds to ', vector_counter
    print 'Left with ', read_number - vector_counter        