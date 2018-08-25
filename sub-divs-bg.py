#!/usr/bin/env python
from __future__ import division
import sys
from collections import defaultdict
import numpy as np

#import timeit

def zero_if_lower_than2(x):
    if x < 2:
        return 0
    return x

def gothrubedgraph(bedgraph):
    bg=open(bedgraph,'r')
    chrcov=defaultdict(list)
    for i,line in enumerate(bg):
        if i%1000000==0:
            print str(i/1000000)+'M','lines read'
        if not line.startswith('track'):
            (ch,start,end,cov)=line.rstrip().split('\t')
            chrcov[ch].append((int(start),int(end),int(cov)))
            #print chr,start,end,cov
    return chrcov

def getchrsize(sfile):
    sizes=defaultdict(int)
    for line in open(sfile):
        (ch,size)=line.rstrip().split('\t')
        print ch,size
        sizes[ch]=int(size)
    return sizes
    
def fill(size,bedlist):
    chrlist=[0]*size
    for pch in bedlist:
        for pos in range(pch[0],pch[1]):
            chrlist[pos]=pch[2]
    return chrlist

def list2bedgraph(ch,plist):
    ppos=0
    pcov=0
    bglist=[]
    for i,pos in enumerate(plist):
        if not pos==pcov:
            bglist.append('\t'.join(map(str,[ch,ppos,i,pcov])))
            pcov=pos
            ppos=i
            
    return bglist

print "Loading file",sys.argv[1] 
bed1=gothrubedgraph(sys.argv[1])
print "Loading file",sys.argv[2] 
bed2=gothrubedgraph(sys.argv[2])
print "Loading sizes",sys.argv[1] 
sizes=getchrsize(sys.argv[3])
#sizes=defaultdict(int)
#sizes['chr1']=18000
print len(sizes.keys()),'chromosomes to deal with...'
#print len(bed1),bed1
out=open(sys.argv[4],'w')
out.write("track type=bedGraph\n")
ct=0
for ch in sizes:
    #print ch
    chr1=fill(sizes[ch],bed1[ch])
    chr2=fill(sizes[ch],bed2[ch])
    #chrd = np.array(chr1)/ np.array(chr2)
    chrd = [c1 / c2 if c2!= 0 else 0 for c1,c2 in zip(chr1, chr2)]
    np_chrd = np.array(chrd)
    max_chrd = max(chrd)
    chrd_cutoff = [zero_if_lower_than2(x) for x in np_chrd]
    chr_sub =list(np.array(chr1) - np.array(chr2))
    new = [c_div if cd_sub > 0 else 0 for c_div,cd_sub in zip(chrd_cutoff, chr_sub)]
    #print chr1[10:100]
    #print chr2[10:100]
    #print chrs[10:100]
    bgs=list2bedgraph(ch,new)
    if len(bgs)>0:
        out.write('\n'.join(bgs)+'\n')
    ct+=1
    if not '_' in ch:
        print ch,'just processed'
    if ct%1000==0:
        print ct,'chromosomes processes'



    
    
    
