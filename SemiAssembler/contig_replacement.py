#!/usr/bin/env python2.7.9

import sys; 

import pysam


def get_chrom(fasta_fh, chrom):
    return fasta_fh.fetch(chrom)

infile = pysam.AlignmentFile("contig_sort.bam", "rb")

#outfile = pysam.AlignmentFile("head_out", "w", template=infile)
target = open("Target_genome.fa", 'w')
fasta = pysam.Fastafile("target_s1")

target_I = open("Long_inertion", 'w')
target_D = open("Long_deletion", 'w')

insertion_count= open("insertion_count", 'w')
deletion_count= open("deletion_count", 'w')

unmapped_len_count= open("unmapped_len", 'w')

#infile = pysam.AlignmentFile("t.bam", "rb")
#target = open("t_python", 'w')
#fasta = pysam.Fastafile("t_ref.fa")



flag=1;
i=0;
for read in infile:
    i+=1;
    if(not(read.is_unmapped)):
        cigarLine=read.cigar;
        qend_pos=read.query_alignment_end;
        #print qend_pos;
        #print cigarLine;
        aln_len=0;
        j=0;
        s_count_s=0;
        e_count_s=0;
        cigarLine_count=len(cigarLine) 
        #print cigarLine_count;
        longest_I=0;
        longest_D=0;
        for(cigarType, cigarLength) in cigarLine:
            j+=1;
            if (  cigarType == 0):
                aln_len+=cigarLength;
            elif (  cigarType == 1):
                insertion_count.write(str(cigarLength));
                insertion_count.write("\n");
                longest_I+=cigarLength;
            elif (  cigarType == 2):
                deletion_count.write(str(cigarLength));
                deletion_count.write("\n");
                aln_len+=cigarLength;
                longest_D+=cigarLength;
            elif (  cigarType == 3):
                aln_len+=cigarLength;
            elif (  cigarType == 4):
                    if( j == 1):
                        s_count_s=cigarLength;
                    if( j == cigarLine_count):
                        e_count_s=cigarLength;
        if(longest_I > 0):
            target_I.write(str(longest_I));
            target_I.write("\n");
        if(longest_D > 0):
            target_D.write(str(longest_D));
            target_D.write("\n");
       
        start_pos=read.reference_start+1;
        end_pos=start_pos+aln_len-1;
        flag_1=start_pos-1;
        #print "s_pos=%d, end_pos=%d, flag_1=%d" % (start_pos,end_pos,flag_1);
        #q_seq=read.query_alignment_sequence;
        q_len=len(read.query_sequence) 
        q_seq_my=read.seq[s_count_s:q_len-e_count_s]
        #print "aln_len=%d, sxs=%d, exs=%d seq=%s" % (aln_len,s_count_s,e_count_s,q_seq);
        #print "aln_len=%d, sxs=%d, exs=%d seq=%s" % (aln_len,s_count_s,e_count_s,q_seq_my);
        if(i==1):
            chr_pre=read.reference_id;
            curr_chrom_name = infile.getrname(read.tid)
            curr_chrom_length= fasta.get_reference_length(curr_chrom_name)
           
            #print curr_chrom_name;
            target.write(">")
            target.write(curr_chrom_name)
            target.write("\n")
            r_seq = fasta.fetch(curr_chrom_name, 0, flag_1)
            target.write(r_seq)
            target.write(q_seq_my)
            flag=end_pos
            #print "r_seq=%s, q_seq=%s, flag=%d" % (r_seq,q_seq_my,flag);
            #curr_chrom_seq = get_chrom(fasta, read.getrname(curr_chrom_id))
        if(i>1):
            chr_curr=read.reference_id;
            if(chr_curr != chr_pre):
                #print "pre=%s, curr=%s" % (chr_pre,chr_curr)
                r_seq = fasta.fetch(curr_chrom_name, flag, curr_chrom_length)
                target.write(r_seq)
                curr_chrom_name = infile.getrname(read.tid)
                curr_chrom_length= fasta.get_reference_length(curr_chrom_name)
                chr_pre=chr_curr;
                #print curr_chrom_name;
                target.write("\n")
                target.write(">")
                target.write(curr_chrom_name)
                target.write("\n")
                r_seq = fasta.fetch(curr_chrom_name, 0, flag_1)
                target.write(r_seq)
                target.write(q_seq_my)
                flag=end_pos;
                #print "r_seq=%s, q_seq=%s, flag=%d" % (r_seq,q_seq_my,flag);
            else:
                if((start_pos <= flag) and (end_pos >= flag)):
                    q_len=len(read.query_sequence)
                    q_seq_my=read.seq[flag-start_pos+s_count_s+1:q_len-e_count_s]
                    target.write(q_seq_my)
                    flag=end_pos;
                    #print "s<f e>f q_seq=%s, flag=%d" % (q_seq_my,flag);
                elif((start_pos < flag) and (end_pos < flag)):
                    continue
                else:
                    r_seq = fasta.fetch(curr_chrom_name, flag, flag_1)
                    target.write(r_seq)
                    target.write(q_seq_my)
                    flag=end_pos;
                    #print "s>f r_seq=%s, q_seq=%s, flag=%d" % (r_seq,q_seq_my,flag);
    else:
        unmapped_len=len(read.query_sequence)
        unmapped_len_count.write(str(unmapped_len))
        unmapped_len_count.write("\n")
r_seq = fasta.fetch(curr_chrom_name, flag, curr_chrom_length)
target.write(r_seq)
target.close()
target_I.close()
target_D.close()
insertion_count.close()
deletion_count.close()
unmapped_len_count.close()
