#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - fixed 'while line.startswith("@"): continue' infinte loop
# - for contigs with non-unique reads, add to unmapped not
#   to unmapped_reads in gene_map()

# NOTES:
# - some contigs/reads may match to the DB multiple times.

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
from collections import Counter

DNA_DB= sys.argv[1]             # INPUT: DNA db used for BWA alignement
contig2read_file= sys.argv[2]   # INPUT: [contigID, #reads, readIDs ...]
gene2read_file= sys.argv[3]     # OUTPUT: [BWA-aligned geneID, length, #reads, readIDs ...]

# make initial dict of contigID<->readsID(s):
contig2read_map= {}
contig_reads= []                                    # list of just reads
with open(contig2read_file,"r") as mapping:
    for line in mapping:
        if len(line)>5:                             # line starts with 'NODE_'
            entry= line.strip("\n").split("\t")     # break tab-separated into list
            contig2read_map[entry[0]]= entry[2:]    # key=contigID, value=list of readID(s)
            contig_reads.extend(entry[2:])          # append all the reads

# make new dict only of contigs with unique reads:
# (hard to tell w BWA what contigs match better, so for reads associated
# w multiple matched contigs, avoid choosing btw contigs for now.)
contig_reads_count= Counter(contig_reads)           # dict of read<->no. of contigs
contig2read_map_uniq= {}
for contig in contig2read_map:
    for read in contig2read_map[contig]:            # If contig has
        if contig_reads_count[read]>1:              #  a read assoc. w multiple contigs
            break                                   #  then throw the contig away,
    else:
        contig2read_map_uniq[contig]= contig2read_map[contig]
                                                    #  else, store it in the unique dict.
# tracking BWA-assigned:
gene2read_map= {}           # dict of BWA-aligned geneID<->readID(s)
mapped_reads= set()         # tracks BWA-assigned reads
prev_mapping_count= 0

#####################################
# FUNCTION:
# add BWA-aligned reads that meet threshold
# to the aligned geneID<->readID(s) dict:

# additional filtering steps:
# (1) only use contigs that contain unique reads (not shared w other contigs)
# (2) make sure BWA-aligned (.sam file has aligned & non-aligned data)
# (3) only use contigs/reads where >90% length of the seq was matched

# CIGAR string describes how read aligns with the ref. Consists of >=1 components.
# Each component comprises an operator and no. bases which the op applies to.
#
# Operators:
# D	Deletion; the nucleotide is present in the reference but not in the read.
# H	Hard Clipping; the clipped nucleotides are not present in the read.
# I	Insertion; the nucleotide is present in the read  but not in the reference.
# M	Match; can be either an alignment match or mismatch. The nucleotide is present in the reference.
# N	Skipped region; a region of nucleotides is not present in the read.
# P	Padding; padded area in the read and not in the reference.
# S	Soft Clipping;  the clipped nucleotides are present in the read.
# =	Read Match; the nucleotide is present in the reference.
# X	Read Mismatch; the nucleotide is present in the reference.


# call as gene_map(BWA .sam file, list of unmapped readIDs)
def gene_map(sam, unmapped):

    len_chars= ["M","I","S","=","X"]    # These particular CIGAR operation cause the
                                        #  alignment to step along the query sequence.
                                        #  Sum of lengths of these operations=length of seq.

    # process .sam file, one contig/read (query) at a time:
    with open(sam,"r") as samfile:
        for line in samfile:
        
            # extract & store data:
            if line.startswith("@") or len(line)<=1:    # If length of line <=1 or line is a header (@...)
                continue                                #  go to the next query (restart for).
            line_parts= line.split("\t")                # Otherwise, split into tab-delimited fields and store:
            query= line_parts[0]                        #  queryID= contig/readID,
            db_match= line_parts[2]                     #  geneID, and a
            flag= bin(int(line_parts[1]))[2:].zfill(11) #  flag---after conversion into 11-digit binary format
                                                        #  where each bit is a flag for a specific descriptor.
            # check if contig:
            if query in contig2read_map:                # If query is a contig (searches through keys)
                if query in contig2read_map_uniq:       #  and it's made of contig-unique reads,
                    contig= True                        #  then mark as contig and continue.
                else:                                   # Otherwise, contig contains non-unique reads,
                    unmapped.add(query)                 #  therefore add it to the unmapped set and
                    continue                            #  go to the next query.
            else:
                contig= False                           # If query isn't a contig, just move on.
                
            # is contig/read BWA aligned?
            if flag[8]=="0":                            # If contig/read is BWA ALIGNED (9th digit=0):
            
                # extract CIGAR data:
                CIGAR= re.split("([MIDNSHPX=])", line_parts[5])
                                                        # Split CIGAR string into list, placing all chars
                                                        #  w/in [...] into own field
                                                        #  (e.g., 9S41M50S->['9','S','41','M','50','S','']).
                # test threshold:
                length= 0
                matched= 0
                for index in range(len(CIGAR))[:-1]:    # Loop CIGAR elements (last element=''),
                    if CIGAR[index+1] in len_chars:     # Use CIGAR operations that step along the query seq,
                        length+= int(CIGAR[index])      #  to determine length of query.
                    if CIGAR[index+1]=="M":             # Use CIGAR match operation to
                        matched+= int(CIGAR[index])     #  determine no. nuclotides matched.
                if matched>length*0.9:                  # If alignment is >90% matched:
                
                    # RECORD alignment:
                    
                    # query aligns to gene with previous alignment(s):
                    if db_match in gene2read_map:                       # If query matches to a previously found gene:
                        if contig:                                      # If query is a contig, then
                            for read in contig2read_map_uniq[query]:    #  take all reads making up that contig
                                gene2read_map[db_match].append(read)    #  and append their readIDs to aligned gene<->read dict,
                                mapped_reads.add(read)                  #  and mark them as assigned by BWA.
                        elif not contig:                                # If query is a read,
                            if query in contig_reads:                   #  but it's part of any another contig,
                                continue                                #  skip to the next query w/out adding to unmapped set.
                            elif query in mapped_reads:                 # If the read has already been assigned by BWA
                                unmapped.add(query)                     #  add it to the unmapped set,
                                for gene in gene2read_map:              #  and
                                    if query in gene2read_map[gene]:    #  remove it from any gene to which it has
                                        if len(gene2read_map[gene])==1: #  been previously aligned to...
                                            del gene2read_map[gene]     #  deleting the entire gene, if it's the
                                        else:                           #  only read mapped to that gene...
                                            gene2read_map[gene].remove(query)
                                        break                           #  (break, as this will only happen once, if any).
                            else:                                       # Otherwise,
                                gene2read_map[db_match].append(query)   #  append its readID to aligned gene<->read dict,
                                mapped_reads.add(query)                 #  and mark it as assigned by BWA.

                    # query matches to BWA-novel gene:
                    else:                                               # If query matches to a novel gene:
                        if contig:                                      # If query is a contig, then
                            gene2read_map[db_match]= list(contig2read_map_uniq[query])
                                                                        #  add the readIDs of all the reads making up
                                                                        #  that contig to the aligned gene<->read dict,
                            for read in contig2read_map_uniq[query]:    #  and mark all those reads
                                mapped_reads.add(read)                  #  as assigned by BWA.
                        elif not contig:                                # If query is a read
                            if query in mapped_reads:                   #  but it has already been assigned by BWA,
                                unmapped.add(query)                     #  add it to the unmapped set,
                                for gene in gene2read_map:              #  and
                                    if query in gene2read_map[gene]:    #  remove it from any gene to which it has
                                        if len(gene2read_map[gene])==1: #  been previously aligned to...
                                            del gene2read_map[gene]     #  deleting the entire gene, if it's the
                                        else:                           #  only read mapped to that gene...
                                            gene2read_map[gene].remove(query)
                                        break                           #  (break, as this will only happen once, if any).
                            else:                                       # Otherwise,
                                gene2read_map[db_match]= [query]        #  add its readID to aligned gene<->read dict
                                mapped_reads.add(query)                 #  and mark it as assigned by BWA.
                
                    continue    # If query was BWA-aligned and met threshold
                                #  (and was therefore processed), skip to next query
                    
            unmapped.add(query) # If query was not BWA-aligned or failed match threshold
                                # put queryID back in unmapped set.

#####################################

# loop over remainder (after argv[4]) in sets of 3:
# (readtype sets: contigs, merged, unmerged1, unmerged2)
# (same .sam file for unmerged1 & unmerged2)
for x in range((len(sys.argv)-4)/3):
    read_file= sys.argv[3*x+4]      # INPUT: all contig/readIDs and seqs (.fasta)
    read_seqs= SeqIO.index(read_file, os.path.splitext(read_file)[1][1:])
                                    # dict of all read SeqRecords: key=contig/readID
                                    #  (second argument specifies filetype, e.g., "fasta")
    BWA_sam_file= sys.argv[3*x+5]   # INPUT: BWA-aligned&unaligned contig/readIDs (.sam)
    output_file= sys.argv[3*x+6]    # OUTPUT: non-BWA-aligned contig/readIDs and seqs (.fasta)

    # "non-BWA-aligned" consists of:
    # (1) contigs containing non-unique reads,
    # (2) contigs/reads not aligned by BWA at all, and
    # (3) BWA-aligned contigs/reads that don't meet match % threshold:
    unmapped_reads= set()           # for set of readIDs
    unmapped_seqs= []               # for list of SeqRecords

    # process BWA-aligned reads:
    gene_map(BWA_sam_file, unmapped_reads)
    #  Novel reads that meet threshold, append to gene2read_map (aligned geneID<->readID(s) dict).
    #  contig/readIDs of contigs/reads that never meet threshold are placed in unmapped_reads.
    
    # WRITE OUTPUT: non-BWA-aligned contig/readIDs
    # and seqs (.fasta):
    for read in unmapped_reads:                     # Put corresponding SeqRecords for unmapped_reads
        unmapped_seqs.append(read_seqs[read])       #  into unmapped_seqs
    with open(output_file,"w") as out:
        SeqIO.write(unmapped_seqs, out, "fasta")    #  and write it to file.

    # print no. aligned reads from current readtype set:
    print str(len(mapped_reads)-prev_mapping_count) + " additional reads were mapped from " + os.path.basename(read_file)
    prev_mapping_count= len(mapped_reads)


# WRITE OUTPUT: write gene2read file of BWA-aligned:
# [BWA-aligned geneID, length, #reads, readIDs ...]
genes= []
with open(gene2read_file,"w") as out_map:
    for record in SeqIO.parse(DNA_DB, "fasta"):         # Loop through SeqRec of all genes in DNA db:
                                                        #  (DNA db is needed to get the sequence.)
        if record.id in gene2read_map:                  #  If DNA db gene is one of the matched genes,
            genes.append(record)                        #  append the SeqRec to genes list (for next file), and
            out_map.write(record.id + "\t" + str(len(record.seq)) + "\t" + str(len(gene2read_map[record.id])))
                                                        #  write [aligned geneID, length, #reads],
            for read in gene2read_map[record.id]:
                out_map.write("\t" + read.strip("\n"))  #  [readIDs ...],
            else:
                out_map.write("\n")                     #  and a new line character.

print "Reads mapped to %d genes." % (len(genes))