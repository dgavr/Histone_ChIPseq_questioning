#!/usr/bin/env python

###############################################################
# Goal: Looks for overlap of H3K27ac histone marks and H3Kme1 #
###############################################################


#new !!!! 

import sys
#from collections import defaultdict

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
    bedrecord = dict()          #Create dictionary with information from this bed file
    #print 'Hi I just made a dictionary in which to parse the contents of your bed file'
    beddata = open(bed_file,'r')
    #print 'OK, so now I am opening your  file called', bed_file
    thereadline= beddata.readline()
    #print thereadline
    if "track" in thereadline : 
        #if len(peak.split())>3
        for peak in beddata: 
            if peak.split()[0] in bedrecord:        #if record is not already present in dictionary
                #print 'Another peak on it is located at', peak.split()[0] 
                bedrecord[peak.split()[0]].append("-".join(peak.split()[1:3])) 
            else : 
                #bedrecord[peak.split()[0]]= str(peak.split()[1],'-',peak.split()[2])
                bedrecord[peak.split()[0]]= ["-".join(peak.split()[1:3])]           #if it IS already present, just add onto existing record
                #print "New chromosome in your file - it is ", peak.split()[0] 
    else:
        print 'This is not a bed file dude!'
    return bedrecord


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

def convert_dict_to_bed(dicti, chr_list,file_name):
    """
    (dict) --> (file)
    Writes contents of a dictionary to file 
    """
    filename= file_name+'.bed'
    myFile = open(filename,'w')
    myFile.write('track name=Inter'+sys.argv[1]+'-'+sys.argv[2]+'\n')
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


# def positive_complicated_fxn(track1_list_element, track2_list_element):
#     """ 
#     (list), (list) --> (string)
#     """
#     start_track1 = int(track1_list_element.split('-')[0])
#     end_track1 = int(track1_list_element.split('-')[1])
#     start_track2 = int(track2_list_element.split('-')[0])
#     end_track2 = int(track2_list_element.split('-')[1])
# 
#     if not ((start_track1 >= end_track2) or (start_track2 >= end_track1 )):   
#         max_start = max(start_track1, start_track2)
#         min_end = min(end_track1, end_track2)
#         result = str(max_start)+'-'+str(min_end)
#         #print track1_list_element,track2_list_element," : ",result
#         return result


def overlap(list1,list2):
    """"
    (list)(list)->(list)(list)
    Takes two lists consisting of start and end positions in bed file for each histone mark track,
    and concatenates them into a single list, then sorts them first by start or end and then by value of position. 
    """
    
    coord=[]
    for pos1 in list1:
        #print 'pos in list1 is', pos1
        coord.append(('S',int(pos1.split('-')[0]), 'l1'))
        #print 'S is ', pos1.split('-')[0]
        coord.append(('E',int(pos1.split('-')[1]),'l1'))
        #print 'E is ', pos1.split('-')[1]
        #print coord 
    for pos2 in list2:
        #print 'pos in list2 is', pos2
        coord.append(('S',int(pos2.split('-')[0]),'l2'))
        #print 'S is ', pos2.split('-')[0]
        coord.append(('E', int(pos2.split('-')[1]),'l2'))
        #print 'E is ', pos2.split('-')[1]
        #print coord
        
    coord.sort(key = lambda x : x[0], reverse = True)
    #print 'coord after first sort \n', coord
    coord.sort(key = lambda x : x[1])
    #print 'coord after 2nd sort by number \n', coord
    # PART 1: SEARCHES FOR OVERLAPS BETWEEN 2 HISTONE MARKS
    new_coord_list = []     #initialize new list to which to move all those that don't overlap
    #index = 0       #position in list 
    spos=0          # start pos initialized 
    ct=0
    ovl=[]
    for pos in coord:
        new_coord_list.append(pos)
        #print pos, 'doesn\'t overlap'
        index = int(new_coord_list.index(pos)) 
        if pos[0]=='S':
            ct+=1
            if ct==2:
                spos=pos[1]
        if pos[0]=='E':
            ct-=1
            if ct==1:
                if not spos==pos[1]:
                    #print spos, '-', pos[1], 'overlap'
                    ovl.append(('ovl', spos, pos[1])) # add to overlap vector the positions that overlap
                    #print 'overlap found! :', [str(spos),str(pos[1]),'ovl']
                    #print 'removing ', new_coord_list[index]
                    del new_coord_list[index]
                    #print 'removing', new_coord_list[index-1]
                    del new_coord_list[index-1]
                    
    #                
    new_coord_list.sort(key = lambda x : x[0], reverse = True)
    start=0
    end = 0
    two_hist_away_from_cent_of_peak = 0
    two_hist_away_list = []
    for nc_pos in new_coord_list:
        if nc_pos[0]=='S':
            if (start<=two_hist_away_from_cent_of_peak) and (two_hist_away_from_cent_of_peak !=0) and (end!=0): 
            #if center_of_peak <= two_hist_away_from_cent_of_peak and (two_hist_away_from_cent_of_peak !=0):
                two_hist_away_list.append('-'.join([str(start),str(end), 'tha']))
                start= nc_pos[1]
        if nc_pos[0]=='E':
            end = nc_pos[1]
            center_of_peak= (start+nc_pos[1])/2
            two_hist_away_from_cent_of_peak = center_of_peak + 300
    # print 'new_coord_list: ', new_coord_list
    return ovl, new_coord_list
  
def parse_overlaps_and_two_hist_away(ovl_list, new_coord_list):
    """
    (list)(list) -> (string)
    
    """
    both_list = new_coord_list[:]
    # print 'before adding overlaps both list looks like this'
    # print both_list
    for ovl_element in ovl_list: 
        both_list.append(ovl_element)
    # print 'after adding overlaps, both list looks like this'
    # print both_list 
    
    both_list.sort(key = lambda x : x[0], reverse = True)
    both_list.sort(key = lambda x : x[1])
    # print 'after sorting both lists looks like this'
    # print both_list 

    dist_list = []
    final_peak_list = []
    case1_counter = 0
    case2_counter = 0
    case3_counter = 0
    case4_counter = 0
    #print 'size is ', len(both_list)
    #range = range(len(both_list))
    for i in range(len(both_list)):
        #print 'both_list[i] is', both_list[i]
        #########################
        #case1: ovl --500bp -4kb -- ovl #
        #########################
        if i+1 < len(both_list): 
            if 'ovl' in both_list[i][0] and 'ovl' in both_list[i+1][0]:
                #print 'case1 ovl --?--- ovl '
                beg_nfr =both_list[i][2]      #end of first peak
                #print 'beg_nfr is', beg_nfr
                end_nfr =both_list[i+1][1] 
                #print 'end_nfr is', end_nfr    #beginning of second peak
                dist_bet_peaks = int(end_nfr) - int(beg_nfr)
                #print 'nfr is ', dist_bet_peaks
                if dist_bet_peaks > 500 and dist_bet_peaks < 4000: 
                    #print 'yay right size nfr - add to list'
                    final_peak_list.append(str(both_list[i][1])+'-'+str(both_list[i][2]))
                    final_peak_list.append(str(both_list[i+1][1])+'-'+str(both_list[i+1][2]))
                    #print 'so it looks like this:', final_peak_list
                    dist_list.append('-'.join([str(beg_nfr),str(end_nfr)]))
                    case1_counter+=1
                #else:
                    #print 'not the right size'
        #################################       
        #case2: ovl -- --500bp -4kb --  l1/2 --l1/2 #
        #################################
        if i+4 < len(both_list):
            if 'ovl' in both_list[i][0] and ('S' in both_list[i+1][0]) and ('E' in both_list[i+2][0]) and 'S' in both_list[i+3][0] and 'E' in both_list[i+4][0]:
                    #print 'case2 ovl --? --- l1/2 l1/2 '
                    if both_list[i+1][2] != both_list[i+3][2]:      #making sure two consecutive peaks are different histone marks rather than two consecutive histone marks (eg. 2x 27ac) 
                        beg_nfr = both_list[i][2]        #end of first peak
                        end_nfr = both_list[i+1][1]                 #beginning of next peak
                        dist_bet_peaks =  int(end_nfr) - int(beg_nfr)               
                        if dist_bet_peaks > 500 and dist_bet_peaks < 4000:
                            final_peak_list.append(str(both_list[i][1])+'-'+str(both_list[i][2])) # overlap peak
                            final_peak_list.append(str(both_list[i+1][1])+'-'+str(both_list[i+2][1])) #  non-overlapping peak on other side of nfr
                            final_peak_list.append(str(both_list[i+3][1])+'-'+str(both_list[i+4][1])) # second consecutive peak on the other side of 
                            dist_list.append('-'.join([str(beg_nfr),str(end_nfr)]))
                            case2_counter+=1
        ####################################
        #case3: l1/2 --- l1/2 --500bp -4kb --  ovl #
        ####################################
            if ('S' in both_list[i][0])  and ('E' in both_list[i+1][0]) and ('S' in both_list[i+2][0])  and ('E' in both_list[i+3][0]) and  'ovl' in both_list[i+4][0]:
                if both_list[i+1][2] != both_list[i+3][2]:      #making sure two consecutive peaks are different histone marks rather than two consecutive histone marks (eg. 2x 27ac) 
                    beg_nfr = both_list[i+3][1]             #end of first peak
                    end_nfr = both_list[i+4][2]               #beginning of next peak
                    dist_bet_peaks = int(end_nfr) - int(beg_nfr)
                    if dist_bet_peaks > 500 and dist_bet_peaks < 4000:
                        final_peak_list.append(str(both_list[i][1])+'-'+str(both_list[i+1][1]))
                        final_peak_list.append(str(both_list[i+2][1])+'-'+str(both_list[i+3][1]))
                        final_peak_list.append(str(both_list[i+4][1])+'-'+str(both_list[i+4][2]))
                        dist_list.append('-'.join([str(beg_nfr), str(end_nfr)]))
                        case3_counter +=1
        ###########################################        
        #case4: l1/2 -- l1/2 --500bp -4kb --  l1/2 --l1/2 #
        ###########################################
        if i+7 < len(both_list):
            if ('S' in both_list[i][0])  and ('E' in both_list[i+1][0]) and ('S' in both_list[i+2][0])  and ('E' in both_list[i+3][0])and ('S' in both_list[i+4][0]) and ('E' in both_list[i+5][0]) and ('S' in both_list[i+6][0]) and ('E'in both_list[i+7][0]):
                if both_list[i][0] != both_list[i+2][0] and both_list[i+4] != both_list[i+6]:
                    beg_nfr = both_list[i+3][1]
                    end_nfr = both_list[i+4][1]
                    dist_bet_peaks = int(end_nfr) - int(beg_nfr)
                    if dist_bet_peaks > 500 and dist_bet_peaks < 4000:
                        final_peak_list.append(str(both_list[i][1]+'-'+str(both_list[i+1][1]))) # first n-o peak
                        final_peak_list.append(str(both_list[i+2][1]+'-'+str(both_list[i+3][1]))) #2nd n-o peak
                        final_peak_list.append(str(both_list[i+4][1]+'-'+str(both_list[i+5][1]))) #3rd n-o peak
                        final_peak_list.append(str(both_list[i+6][1]+'-'+str(both_list[i+7][1]))) #4th n-o peak
                        dist_list.append('-'.join([str(beg_nfr),str(end_nfr)]))
                        case4_counter +=1
    # print 'Number of case 1 scenarios:', case1_counter
    # print 'Number of case 2 scenarios:', case2_counter
    # print 'Number of case 3 scenarios:', case3_counter
    # print 'Number of case 4 scenarios:', case4_counter
    # print 'dist_list is', dist_list
    return dist_list, final_peak_list


if __name__ == "__main__":
    """ Opens several bed files, creates a dictionary with the overlaps and
         writes dictionary to file in bed format.  
    """
    ftrack1 = sys.argv[1]   #Track file (bed format) input
    ftrack2 = sys.argv[2]   #Track file (bed format) input
    name = ftrack1.split('.')[0]+'_'+ftrack2.split('.')[0] #Name to add to NEW OUTPUT file 
    #print ftrack1,ftrack2
    nucleosome_depleted_dict = dict()  #initialize dictionary
    final_peak_dict = dict()
    
    #trackp1 = go_thru('GSM803827_H3K4me1_80_epi_peaks.bed')
    #trackp2 = go_thru('GSM803831_H3K27ac_80_epi_peaks.bed')
    trackp1 = go_thru(ftrack1)
    trackp2 = go_thru(ftrack2)
    list_shared_chr = compare_chr(trackp1, trackp2)
    #print len(list_shared_chr),'chromosomes to check...'
    #print list_shared_chr
            #print 'The list of shared chromosomes between the two lists is', list_shared_chr

    ########################################
    # 1. Find overlap between me1 and Ac27 #
    ########################################
    # test1=['41-131','4061-4141','5323-6550','7222-7459','8752-9001']
    # test2=['3928-4092','4367-4852','5520-5772','5985-6100','8997-9111']
    
    #print 'Now I am going to go through the different dictionaries made out of your files ....'
    nbchch=0
    for chr in list_shared_chr:            #Loop through the list of shared chromosomes
        #print 'chr no is', chr
        nbchch+=1
        if nbchch%1000==0:
            print nbchch,'chr processed!'     
        ovl, new_coord_list =  overlap(trackp1[chr],trackp2[chr])
        nucleosome_depleted_list, final_peak_list = parse_overlaps_and_two_hist_away(ovl, new_coord_list)
        nucleosome_depleted_dict[chr]= nucleosome_depleted_list
        final_peak_dict[chr] = final_peak_list
        if chr[0]=='c':
            print chr,':',len(trackp1[chr]),'peaks in', ftrack1,',', len(trackp2[chr]),'peaks in ', ftrack2, len(nucleosome_depleted_list),' overlaps!'
    convert_dict_to_bed(nucleosome_depleted_dict, list_shared_chr, 'E1b_ndr_output_'+name)
    convert_dict_to_bed(final_peak_dict, list_shared_chr, 'E1b_fp_output_'+name)