#!/usr/bin/env python
import sys


# Searches for nucleosome depleted regions in subtraction file produced by subtract_input.py
# Returns bed file with beginning and end of nucleosome depleted regions


if __name__ == "__main__":
    subtracted_filename = sys.argv[1]
    output_file= open('ndr_'+subtracted_filename[subtracted_filename.find('sub')+4:]+'.bed','w')
    chrom = subtracted_filename.split('_')[1]
    counter= 0
    prev_me1 = None
    prev_ac27 = None
    ndr_start = None
    ndr_end = None
    for line in open(subtracted_filename):
        if not line.startswith('Chromosome'):
            me1 = int(line.split('\t')[0])
            ac27 = int(line.split('\t')[1])
            print 'me1 = ', me1, ' prev me1 = ', prev_me1
            print 'ac27 = ', ac27, ' prev_ac27 = ', prev_ac27
            # NDR : both values are negative or = 0
            if (type(me1)!=None) and (type(ac27)!=None) and (type(prev_me1)!=None) and (type(prev_ac27)!=None): 
                if (me1 < 1) and (ac27 < 1):
                    print 'NDR!!!'
                # Start of chromosome, Start NDR
                    if ((type(prev_me1) == None) and  (type(prev_ac27) == None)) or ((prev_me1 > 1) and (prev_ac27 > 1))or ((prev_me1 < 1) and (prev_ac27 > 1)) or ((prev_me1 > 1) and (prev_ac27 < 1)):
                        print 'start of NDR '
                        ndr_start = counter
                    elif ((prev_me1 < 1) and (prev_ac27 < 1)) and ((not type(prev_me1)==None) and (not type(prev_ac27) == None)) :
                        print 'Already in NDR'
            if (me1 > 1) and (ac27 > 1) and not ((type(prev_me1)==None) and not (type(prev_ac27) == None)):
                    print 'Not in NDR'
                    if ((prev_me1 < 1) and (prev_ac27 < 1)) and (type(ndr_start)!=None):
                        ndr_end = counter
                        print 'End of NDR!!!!!'
                        print 'chrom is ', type(chrom), 'and is ', chrom
                        print 'ndr_start is', type(ndr_start), 'and is ', ndr_start
                        print 'ndr_end is', type(ndr_end), 'and is', ndr_end 
                        output_file.write('\t'.join([chrom, ndr_start, ndr_end]))
            counter+=1
            prev_me1 = me1
            prev_ac27 = ac27