#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - rewrote contig_map(), also added a check to remove mapped reads from the unmapped set.

# NOTES:
# - For unmerged pairs, where only one of a pair is mapped to a contig, the other
#   will be discarded, and the read considered to be part of the contig.
# - It is possible for a read to align to multiple contigs.

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

# files + extraction:
paired_file1= sys.argv[1]                           # INPUT: unmerged1 readIDs and seqs (.fastq)
paired1_seqs= SeqIO.index(paired_file1, "fastq")    # dict of unmerged1 read SeqRecords: key=readID

paired_file2= sys.argv[2]                           # INPUT: unmerged2 readIDs and seqs (.fastq)
paired2_seqs= SeqIO.index(paired_file2, "fastq")    # dict of unmerged2 read SeqRecords: key=readID

unpaired_file= sys.argv[3]                          # INPUT: merged readIDs and seqs (.fastq)
unpaired_seqs= SeqIO.index(unpaired_file, "fastq")  # dict of merged read SeqRecords: key=readID

paired_sam_file= sys.argv[4]        # INPUT: unmerged reads that are BWA-aligned&unaligned to contigs (.sam)
unpaired_sam_file= sys.argv[5]      # INPUT: merged reads that are BWA-aligned&unaligned to contigs (.sam)

output_file= sys.argv[6]            # OUTPUT: contig<->read map: [contigID, #reads, readIDs ...]

# tracking BWA-aligned:
contig2read_map= defaultdict(set)   # dict of BWA-aligned contigID<->readID(s)
mapped_reads= set()                 # tracks BWA-assigned reads
mapped_list= []

#####################################
# FUNCTION:
def contig_map(sam, unmapped):

    # process .sam file, one contig/read (query) at a time:
    with open(sam,"r") as samfile:
        for line in samfile:
        
            # valid line?
            if line.startswith("@") or len(line)<=1:    # If length of line <=1 or line is a header (@...)
                continue                                #  skip to the next line.
        
            # extract data:
            line_parts= line.strip("\n").split("\t")    # Otherwise, split into tab-delimited fields and extract:
            readID= line_parts[0]                       #  the readID,
            flag= bin(int(line_parts[1]))[2:].zfill(11) #  the flag--converting to 11-digit binary:
                                                        #  each bit is a flag for a specific descriptor.
            # is read BWA-aligned?:
            if flag[8]=="1":                            # If read is NOT BWA-ALIGNED (9th digit=1),
                unmapped.add(readID)                    #  add it to the unmapped set and
                continue                                #  skip to the next read.
        
            # store info in map:
            contigID= line_parts[2]                     # Extract the contigID, and
            contig2read_map[contigID].add(readID)       #  add the readID to the contigID dict entry
                                                        #  (it's a set, so it won't store duplicates, though it
                                                        #  is possible for a read to align to multiple contigs),
            mapped_reads.add(readID)                    #  and mark the read as mapped to a contig.

    # Remove reads previously added to the unmapped set but later found to align to a contig:
    # Occurs for unmerged pairs, where only one of a pair is mapped to a contig, the other will
    # need to be removed from the unmapped set.
    print '\n' + str(os.path.basename(sam))
    print 'umapped no. (before double-checking mapped set)= ' + str(len(unmapped))
    for read in mapped_reads:                           # Take all mapped reads and
        try:                                            #  if the exist in the unmapped set, then
            unmapped.remove(read)                       #  remove them from the unmapped set.
        except:
            pass
    print 'umapped no. (after double-checking mapped set)= ' + str(len(unmapped))

#####################################

# extraction of contig-aligned reads and unaligned reads:
unmapped_paired_reads= set()
unmapped_unpaired_reads= set()
contig_map(paired_sam_file, unmapped_paired_reads)      # unmerged reads
contig_map(unpaired_sam_file, unmapped_unpaired_reads)  # merged reads

# WRITE OUTPUT: non-contig unmerged:
unmapped_paired1_seqs= []
unmapped_paired2_seqs= []
for read in unmapped_paired_reads:                      # For the non-contig unmerged reads,
    unmapped_paired1_seqs.append(paired1_seqs[read])    #  isolate the corresponding SeqRecord for fwd read,
    unmapped_paired2_seqs.append(paired2_seqs[read])    #  and for the rev read.
with open(os.path.splitext(paired_file1)[0] + "_unmapped.fastq","w") as out:
    SeqIO.write(unmapped_paired1_seqs, out, "fastq")    # Write the sequence for fwd read to file
with open(os.path.splitext(paired_file2)[0] + "_unmapped.fastq","w") as out:
    SeqIO.write(unmapped_paired2_seqs, out, "fastq")    #  and for the rev read.

# WRITE OUTPUT: non-contig merged
unmapped_unpaired_seqs= []
for read in unmapped_unpaired_reads:                    # For the non-contig merged reads,
    unmapped_unpaired_seqs.append(unpaired_seqs[read])  #  isolate the corresponding SeqRecord,
with open(os.path.splitext(unpaired_file)[0] + "_unmapped.fastq","w") as out:
    SeqIO.write(unmapped_unpaired_seqs, out, "fastq")   #  and write the sequence to file.

# WRITE OUTPUT: contig<->read map
# [contigID, #reads in contig, read IDs...]
with open(output_file,"w") as outfile:
    for contig in contig2read_map:                      # For each contig,
        outfile.write(contig + "\t" + str(len(contig2read_map[contig])))
                                                        # write [contigID, #reads in contig, ...]
        for read in list(contig2read_map[contig]):
            outfile.write("\t" + read)                  # [readIDs ...],
        else:
            outfile.write("\n")                         # and a new line character.
