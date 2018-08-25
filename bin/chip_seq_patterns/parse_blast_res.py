#!/usr/bin/env python

import sys
from collections import defaultdict


def make_header_dict(fasta_file):
    """
    """
    header_list_dict = defaultdict(list)
    for line in open(fasta_file):
        if line.startswith('>'):
            header_list_dict[line[1:].split(' ')[0]]= line[1:].split(' ')[1]
    return header_list_dict

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    blast_file = sys.argv[2]
    out = open('names_'+blast_file, 'w')
    fa_header_dict = make_header_dict(fasta_file)
    for line in open(blast_file):
        (line) = line.split('\t')
        print line[0]
        if line[0] in fa_header_dict: 
            name = fa_header_dict[line[0]]
        print name
        out.write(str(name)+'\t'+ '\t'.join([str(line[1:])])+'\n')
        if line[0] not in fa_header_dict:
            print str(line[0]), ' - not in dictionary!'
            