#!/usr/bin/env python2.7.9

import sys; 
import os;
import time
import pysam
import numpy as np

fname_pe="pe_sort.bam"
fname_mt="mt_sort.bam"
read_coverage=50

thread_t=read_coverage*0.08+1
xsxm_count=0;
xsxm_thread=thread_t;
xmxs_thread=thread_t;


infile = pysam.AlignmentFile(fname_pe, "rb")
in_mate = pysam.AlignmentFile(fname_mt,"rb")


xsxm_list= open("xsxm_list", 'w')
mate_list= open("mate_list", 'w')
mate_list_D= open("mate_list_D", 'w')
time_list= open("time_count",'w')
find_I_list=open("find_list_I",'w')
find_D_list=open("find_list_D",'w')
###########################################################

def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst)+1)/2)-1]
    else:
            return int(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0


# count average insert size
t_start=time.time()
i=0;
sum=0;
insertsize_list = []
for read in in_mate:
    if(not(read.is_unmapped)):
        if((read.tlen > 0) and (read.tlen < 30000)):
            insertsize_list.append(read.tlen)

average = np.average(insertsize_list)
stdev = np.std(insertsize_list)
thread_mate = average-(2)*stdev
thread_mate_D = average+(2)*stdev
#median_size = median(insertsize_list)

print average
print stdev
print thread_mate
#print median_size

t_end=time.time()
total_time= t_end - t_start;
time_list.write(str(total_time));
time_list.write("\n");
time_list.write("\n");

in_mate.close()

##########################################################
# XSXM > thread and determine insertion or deletion
i=0;
t_start=time.time();
for read in infile:
    i+=1;
    if(not(read.is_unmapped)):
        curr_chrom = read.reference_id;
        curr_pos = read.reference_start;
        if(i==1):
            pre_chrom = read.reference_id;
            pre_pos=0;
        if(curr_pos != pre_pos):
            # if XSXM > thread means possible breakpoint ,and then determine insertion or deletion by median insert size of mate pair that across this breakpoint
            if(xsxm_count > xsxm_thread):
                mt_start_pos=pre_pos-15000
                mt_end_pos=pre_pos + 1             
                if(mt_start_pos < 0):
                    mt_start_pos=0;
   
                chr=pre_chrom_name
                in_mate = pysam.AlignmentFile(fname_mt,"rb")
                iter_mt = in_mate.fetch(chr, mt_start_pos, mt_end_pos)
                flag_I=0
                flag_D=0
                mt_insertsize_list = []
                for mtread in iter_mt:
                    if(not(mtread.is_unmapped)):
                        if((mtread.tlen > 0) and (mtread.tlen < 15000) and (mtread.reference_start <= mt_end_pos) and (mtread.next_reference_start >= mt_end_pos)):
                            mt_insertsize_list.append(mtread.tlen)
                        if((abs(mtread.tlen) < thread_mate) and (mtread.tlen != 0)):
                            flag_I+=1;
                        if((abs(mtread.tlen) > thread_mate_D) and (mtread.tlen != 0)):
                            flag_D+=1;
                median_size = np.median(mt_insertsize_list)
                
                #print median_size
                del mt_insertsize_list[:]
                # possible insertion 
                if (median_size < average): 
                     pre_len= average - median_size
                     find_I_list.write(pre_chrom_name);
                     find_I_list.write("\t");
                     find_I_list.write(str(pre_pos+1));
                     find_I_list.write("\t");
                     find_I_list.write(str(pre_len));
                     find_I_list.write("\n");
                # possible deletion
                if (median_size >= average):
                     pre_len= median_size - average
                     find_D_list.write(pre_chrom_name);
                     find_D_list.write("\t");
                     find_D_list.write(str(pre_pos+1));
                     find_D_list.write("\t");
                     find_D_list.write(str(pre_len));
                     find_D_list.write("\n");
                xsxm_count=0;
            else:
                xsxm_count=0;
                cigarLine=read.cigar;
                cigarLine_count=len(cigarLine)
                j=0;
                flag_xs=0;
                for(cigarType, cigarLength) in cigarLine:
                    j+=1;
                    if((cigarType == 4) and (j == 1)):
                        flag_xs=1;
                    if((cigarType == 0) and (j == cigarLine_count) and (flag_xs == 1)):
                        xsxm_count+=1;
                        pre_chrom = read.reference_id;
                        pre_pos = read.reference_start;
                        pre_chrom_name=infile.getrname(read.tid);
        else:
            cigarLine=read.cigar;
            cigarLine_count=len(cigarLine)
            j=0;
            flag_xs=0;
            for(cigarType, cigarLength) in cigarLine:
                j+=1;
                if((cigarType == 4) and (j == 1)):
                    flag_xs=1;
                if((cigarType == 0) and (j == cigarLine_count) and (flag_xs == 1)):
                    xsxm_count+=1;
                    pre_chrom = read.reference_id;
                    pre_pos = read.reference_start;
                    pre_chrom_name=infile.getrname(read.tid);

xsxm_list.close()
infile.close()

t_end=time.time()
total_time= t_end - t_start;
time_list.write(str(total_time));
time_list.write("\n");

##########################################################################################


time_list.close()
mate_list.close()
mate_list_D.close()
time_list.close()
find_I_list.close()
find_D_list.close()



##########################################################################################
# output insertion info and LR_read_info(Left and Right sequence of possible insertions and predict insertion length )

t_start=time.time();
time_list= open("time_findinfo_new",'w')

infile = pysam.AlignmentFile(fname_pe, "rb")
in_mate = pysam.AlignmentFile(fname_mt,"rb")



mt_avg=average
mt_stdev=stdev

insertion_info= open("insertion_info_new", 'w')
LR_read_info= open("LR_read_info", 'w')

f=file("find_list_I",'r')
for line in f:
    items=line.strip().split('\t')
    chr=items[0]
    pos=items[1]
    pe_start_pos=int(pos)-200
    if(pe_start_pos < 0):
        pe_start_pos=0;
    pe_end_pos=int(pos)+100
    pe_pos=int(pos)-1
    iter_pe = infile.fetch(chr, pe_start_pos, pe_end_pos)
    
    mt_start_pos=int(pos)-15000
    if(mt_start_pos < 0):
        mt_start_pos=0;
    mt_end_pos=int(pos)+1
    mt_pos=int(pos)-1
    iter_mt = in_mate.fetch(chr, mt_start_pos, mt_end_pos)


    xsxm_count = 0 
    xmxs_count = 0
    xm_count = 0
    for read in iter_pe:
        if(not(read.is_unmapped)):
            cigarLine=read.cigar;
            cigarLine_count=len(cigarLine)
            j=0;
            if(read.reference_start == pe_pos):
                flag_xs=0;
                for(cigarType, cigarLength) in cigarLine:
                    j+=1;
                    if((cigarType == 4) and (j == 1)):
                        flag_xs=1;
                    if((cigarType == 0) and (j == cigarLine_count) and (flag_xs == 1)):
                        xsxm_count+=1;
            if((read.reference_start > (pe_pos-100)) and (read.reference_start < pe_pos)):
                flag_xm=0;
                for(cigarType, cigarLength) in cigarLine:
                    j+=1;
                    if((cigarType == 0) and (cigarLine_count == 1)):
                       xm_count+=1;
                    if((cigarType == 0) and (j == 1)):
                       flag_xm=1;
                    if((cigarType == 4) and (j == cigarLine_count) and (flag_xm == 1)):
                       xmxs_count+=1;
            if((read.reference_start < (pe_pos-100))):
                pair_flag=0;
                for(cigarType, cigarLength) in cigarLine:
                    if((cigarType == 0) and (cigarLine_count == 1)):
                        pair_flag=1;
                        L_pos=read.reference_start + 1;
                        L_read=read.seq
            if((read.reference_start > pe_pos) and (pair_flag==1)):
                for(cigarType, cigarLength) in cigarLine:
                    if((cigarType == 0) and (cigarLine_count == 1)):
                        pair_flag=2;
                        R_pos=read.reference_start + 1 + 101;
                        R_read=read.seq            
           

  #  print xsxm_count;
  #  print xmxs_count;
  #  print xm_count;
    insertsize_list = [];
    for mt_read in iter_mt:
        if(not(mt_read.is_unmapped)):
            if((mt_read.tlen > 0) and (mt_read.tlen < 15000) and (mt_read.reference_start <= mt_pos) and (mt_read.next_reference_start >= mt_pos)):
            #    print mt_read.tlen
                insertsize_list.append(mt_read.tlen)
    average = np.average(insertsize_list)
    median_size = np.median(insertsize_list)
    insertion_length = mt_avg - median_size;
#    insertion_length = mt_avg - average;
   # print insertion_length;
    
    curr_chrom = read.reference_id;
    chrom_name=infile.getrname(read.tid);
   
    if((xmxs_count > xmxs_thread) and (xsxm_count > xmxs_thread) and (pair_flag == 2) and (xm_count < 2)):
        insertion_info.write(chrom_name+'\t'+pos+'\t'+str(xmxs_count)+'\t'+str(xsxm_count)+'\t'+str(xm_count)+'\t'+str(insertion_length)+'\n');
        LR_read_info.write('>'+chrom_name+'_'+pos+'_'+str(L_pos)+'_'+str(R_pos)+'\n'+L_read+'\n');
        LR_read_info.write('>'+chrom_name+'_'+pos+'_'+str(L_pos)+'_'+str(R_pos)+'\n'+R_read+'\n');
        LR_read_info.write(str(insertion_length)+'\n');
    pair_flag=0;

t_end=time.time()
total_time= t_end - t_start;
time_list.write(str(total_time));
time_list.write("\n");    
      
f.close()
insertion_info.close()
LR_read_info.close()
time_list.close()
infile.close()
in_mate.close()

####################################################################################################################################################
#output deletion range : output_deletion_range

t_start=time.time();
time_list= open("time_findinfo_deletion",'w')

infile = pysam.AlignmentFile(fname_pe, "rb")


deletion_info= open("deletion_info_new", 'w')
deletion_out = open("output_deletion_range", 'w')

#f=file("find_list_D",'r')
f=file("find_list_D",'r')
for line in f:
    items=line.strip().split('\t')
    chr=items[0]
    pos=items[1]
    pe_start_pos=int(pos)-5000
    if(pe_start_pos < 0):
        pe_start_pos=0;
    pe_end_pos=int(pos)+100
    pe_pos=int(pos)-1
    iter_pe = infile.fetch(chr, pe_start_pos, pe_end_pos)
    iter_pe2 = infile.fetch(chr, pe_start_pos, pe_end_pos)    
   # mt_start_pos=int(pos)-15000
   # if(mt_start_pos < 0):
   #     mt_start_pos=0;
   # mt_end_pos=int(pos)+1
   # mt_pos=int(pos)-1
   # iter_mt = in_mate.fetch(chr, mt_start_pos, mt_end_pos)
    
    i=0;
    for read in iter_pe:
        i+=1;
        if(not(read.is_unmapped)):
            cigarLine=read.cigar;
            cigarLine_count=len(cigarLine)
            j=0;
            if(i==1):
                pre_cigarLine=read.cigar;
                pre_cigarLine_count=len(cigarLine)
                pre_pos=read.reference_start
            if(read.reference_start == pe_pos):
                flag_xm=0;
                flag_is_xmxs=0;
                for(cigarType, cigarLength) in pre_cigarLine:
                    j+=1;
                    if((cigarType == 0) and (j == 1)):
                       flag_xm=1;
                       count_pos_m = cigarLength;
                    if((cigarType == 4) and (j == pre_cigarLine_count) and (flag_xm == 1)):
                       flag_is_xmxs=1;
                       delete_start_pos= count_pos_m + pre_pos;
                if(flag_is_xmxs==0):
                    delete_start_pos=pre_pos
                break
            pre_cigarLine=read.cigar;
            pre_cigarLine_count=len(cigarLine)        
            pre_pos=read.reference_start
                
    i=0
    xsxm_count = 0 
    xmxs_count = 0
    xm_count = 0
    xmxs_start=delete_start_pos-100
  
    for read2 in iter_pe2:
        if(not(read2.is_unmapped)):
            cigarLine=read2.cigar;
            cigarLine_count=len(cigarLine)
            j=0;
            if(read2.reference_start == pe_pos):
                flag_xs=0;
                for(cigarType, cigarLength) in cigarLine:
                    j+=1;
                    if((cigarType == 4) and (j == 1)):
                        flag_xs=1;
                    if((cigarType == 0) and (j == cigarLine_count) and (flag_xs == 1)):
                        xsxm_count+=1;
            if((read2.reference_start > (delete_start_pos-100)) and (read2.reference_start < pe_pos)):
                flag_xm=0;
                for(cigarType, cigarLength) in cigarLine:
                    j+=1;
                    if((cigarType == 0) and (cigarLine_count == 1)):
                       xm_count+=1;
                    if((cigarType == 0) and (j == 1)):
                       flag_xm=1;
                    if((cigarType == 4) and (j == cigarLine_count) and (flag_xm == 1)):
                       xmxs_count+=1;
  
    real_deletion_len=int(pos)-delete_start_pos-1
   # print insertion_length;
    
    curr_chrom = read.reference_id;
    chrom_name=infile.getrname(read.tid);
   # info_list=[xmxs_count,xsxm_count,xm_count,insertion_length]
  
    if((xmxs_count > xmxs_thread) and (xsxm_count > xmxs_thread) and (xm_count < 2) and (real_deletion_len > 0)):
        deletion_info.write(chrom_name+'\t'+pos+'\t'+str(xmxs_count)+'\t'+str(xsxm_count)+'\t'+str(xm_count)+'\t'+str(real_deletion_len)+'\t'+str(delete_start_pos+1)+'-'+pos+'\n');
        deletion_out.write(chrom_name+'\t'+str(delete_start_pos+1)+'\t'+pos+'\tD\tX'+'\n'); 
t_end=time.time()
total_time= t_end - t_start;
time_list.write(str(total_time));
time_list.write("\n");    
      
f.close()
deletion_info.close()
