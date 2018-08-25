#!/usr/bin/env python
#Goes through Me3 bed and looks closest GP in proximity 

import sys


def go_thru(bed_file):
    """(bed file) -> (dictionary)
    Goes through a bed file and generates a dictionary (bed record) containing the
    as the key - chromosome name and as the value - a list of strings each containing
    the start and end values. 
    Example of contents of bed file:
    # track type="bed" name="H3K4me3_pk_80\%epi" converted_by="track.serialize"
    Zv8_NA9    676    4438
    Zv8_NA15    1196    1680
    Zv8_NA16    2767    3031
    Zv8_NA16    3631    5050
    Zv8_NA23    0    29
    Dictionary key = 'Zv8_NA9' 
    Dictionary values = ['676-4438',1196-1680 ... ] 
    """
    bed_dict = {}          #Create dictionary with information from this bed file
    beddata = open(bed_file,'r')
    thereadline= beddata.readline()
    if "track" in thereadline : 
        for peak in beddata: 
            if peak.split()[0] in bed_dict:        #if record is not already present in dictionary
                bed_dict[peak.split()[0]].append("-".join(peak.split()[1:3])) 
            else : 
                bed_dict[peak.split()[0]]= ["-".join(peak.split()[1:3])]           
                #if it IS already present, just add onto existing record
    else:
        print 'This is not a bed file dude!'
    return bed_dict

def gene_pred_to_dict():
    """
    gene PredExt format
        string name;        	"Name of gene (usually transcript_id from GTF)"
        string chrom;       	"Chromosome name"
        char[1] strand;     	"+ or - for strand"
        uint txStart;       	"Transcription start position"
        uint txEnd;         	"Transcription end position"
        uint cdsStart;      	"Coding region start"
        uint cdsEnd;        	"Coding region end"
        uint exonCount;     	"Number of exons"
        uint[exonCount] exonStarts; "Exon start positions"
        uint[exonCount] exonEnds;   "Exon end positions"
        uint id;            	"Unique identifier"
        string name2;       	"Alternate name (e.g. gene_id from GTF)"
        string cdsStartStat; 	"enum('none','unk','incmpl','cmpl')"
        string cdsEndStat;   	"enum('none','unk','incmpl','cmpl')"
        lstring exonFrames; 	"Exon frame offsets {0,1,2}"
    
    """
    gene_pred_handler= open("ensGene.txt", 'r')
    gene_pred_dict= {}
    for line in gene_pred_handler:
        #line.split()[0]
        Ensdart_id = line.split()[1] #Ensembl transcript id/ Gene name?
        chrom = line.split()[2] #chromosome id/name
#        strand = line.split()[3]    #strand + or - 
        txStart = line.split()[4]   #Transcription start position
        txEnd = line.split()[5] # Transcription end position
#        cdsStart= line.split()[6]   #Coding region start
#        cdsEnd = line.split()[7]    #Coding region end
#        exonCount = line.split()[8] #Number of exons
#        exonStarts= line.split()[9] #Exon start positions
#        exonEnds = line.split()[10] #Exon end positions
        #dunno = line.split()[11]
#        uid =  str(line.split()[12])#Unique identifier
#        name2 =  line.split()[13]# Alternate name
#        cdsStartStat= line.split()[14] #enum('none','unk','incmpl','cmpl')
#        cdsEndStat = line.split()[15]  #enum('none','unk','incmpl','cmpl') 
#        exonFrames =line.split()[16]    #Exon frame offsets {0,1,2}
        if chrom in gene_pred_dict: 
            gene_pred_dict[chrom].append(str(txStart)+'-'+str(txEnd)+'-'+str(Ensdart_id))
        else: 
            gene_pred_dict[chrom] = []
            gene_pred_dict[chrom]= [txStart+'-'+txEnd+'-'+Ensdart_id]
    return gene_pred_dict

def compare_chr(track1, track2):
    """ Make chromosomes (the keys to the dictionary generated in go_thru() and export
        them as a list. Then compare the list of chromosomes to make sure that they are 
        all present in both bed files. Export the list of K present in both files. Gets 
        rid of chromosomes absent in one of the files.
    """
    intersect = set(track1.keys()).intersection(set(track2.keys()))
    chr_list = list(intersect)
    #print 'comparing chromosome lists'
    return chr_list

def overlap(list1,list2):
    """
    List1 = me3
    List2 = Gene prediction
    """
    coord=[]
    for pos in list1:
        coord.append(('S',int(pos.split('-')[0]),'me3'))
        coord.append(('E',int(pos.split('-')[1]), 'me3'))
    for pos in list2:
        coord.append(('S',int(pos.split('-')[0]),str(pos.split('-')[2])))
        coord.append(('E',int(pos.split('-')[1]),str(pos.split('-')[2])))
    coord.sort(key = lambda x : x[0], reverse = True)
    coord.sort(key = lambda x : x[1])

    spos=0
    ct=0
    ovl=[]
    for pos in coord:
        if pos[0]=='S':
            ct+=1
            if ct==2:
                spos=pos[1]
                gi_spos= pos[2]
        if pos[0]=='E':
            ct-=1
            if ct==1:
                if not spos==pos[1]:
                    e_id = pos[2]
                    gene_id, ok  = parse_gene_id(gi_spos, e_id)
                    if ok ==0:
                        ovl.append(str(spos)+'-'+str(pos[1])+'-'+gene_id)
    return ovl

def parse_gene_id(s_id, e_id):
    """"
    """
    if s_id == 'me3' and e_id!= 'me3':
        gene_id = e_id
        ok = 0
    elif e_id == 'me3' and s_id != 'me3': 
        gene_id = s_id
        ok = 0 
    elif s_id== 'me3' and e_id=='me3':
        gene_id = '2 me3 overlaps?!'
        ok =1
    else:
        gene_id = 'Unknown'
        ok = 1
    return gene_id, ok
    
def convert_dict_to_bed(dicti, chr_list,file_name):
    """
    (dict) --> (file)
    Writes contents of a dictionary to file 
    """
    myFile = open(file_name,'w')
    myFile.write('track name=Me3_expressed'+file_name+'\n')
    #print ' I just made a file to put your results in called'+ filename
    for chr in chr_list: 
        #print 'the chr is', chr
        #print "chr no.",chr," : ",len(dicti[chr])," put enhancers"
        for item in dicti[chr]:
            #print 'the item is ', item
            line=item.split('-')
            line.insert(0,chr)
            #print line
            myFile.write('\t'.join(line)+'\n')
    myFile.close()


if __name__ == "__main__":
    me3_bed_file = sys.argv[1]
    output_file = me3_bed_file.split(".")[0]+'_gpred_output.bed'
    me3_dict = go_thru(me3_bed_file)
    g_pred_dict = gene_pred_to_dict()
    overlap_track ={}
    common_chr = compare_chr(me3_dict, g_pred_dict)
    nbchch=0
    for chr in common_chr: 
        #print 'THIS IS CHROMOSOME NO:', chr
        nbchch+=1
        if nbchch%1000==0:
            print nbchch,'chr processed!'   
        overlap_track[chr] = overlap(me3_dict[chr], g_pred_dict[chr])
    convert_dict_to_bed(overlap_track, common_chr, output_file)