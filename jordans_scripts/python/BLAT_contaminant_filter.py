#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - closed the BLAT_tab_file, after contaminated_seqs is made

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

# INPUTS:
input_file = sys.argv[1]            # INPUT: non-BWA-aligned reads (.fastq); input to BLAT
input_seqs = SeqIO.to_dict(SeqIO.parse(input_file, "fastq")) # dict of non-BWA-aligned read SeqRecords
                                                             #  SeqIO.parse parses the .fastq one read at a time, 
                                                             #  creates a SeqRecord object > SeqIO.to_dict places 
                                                             #  SeqRecord into dictionary using readID as dict key.
BLAT_tab_file = sys.argv[2]         # INPUT: BLAT output (.blatout) (BLAT uses non-BWA aligned reads as input)
                                    #  ****IMPORTANT**** .blatout has identical readIDs on ajacent lines.

# OUTPUTS:
output_file = sys.argv[3]           # OUTPUT: non-BWA&BLAT-aligned reads (.fastq)
output_seqs = []                    # list of non-BWA&BLAT-aligned read SeqRecords
output_file_made = False
contaminat_output_file = sys.argv[4]# OUTPUT: BLAT-aligned reads
contaminat_output_seqs = []         # list of BLAT-aligned read SeqRecords (.fastq)
contaminat_output_file_made = False

# make a list of BLAT-aligned readIDs:
with open(BLAT_tab_file, "r") as tabfile:       # Open .blatout file.
    contaminated_seqs = []
    query_seq = ""
    for line in tabfile:                        # In the .blatout file:
        if len(line) < 2:                       # If length of line < 2,
            continue                            #  go to next line (restart for).
        line_parts = line.split("\t")           # Split line by tab separator (Element 1=readID).
        if query_seq == line_parts[0]:          # If readID is empty or same as last line,
            continue                            #  skip to next line.
        else:
            query_seq = line_parts[0]           # If readID exists, append to a
            contaminated_seqs.append(query_seq) #  list of BLAT-aligned readIDs. 

# separate non-BWA-aliged reads
# according to BLAT-alignment:
for seq in input_seqs:                      # Look through all non-BWA-aligned readIDs:

    # non-BLAT-aligned:                      
    if seq not in contaminated_seqs:        # If read is non-BLAT-aligned,  
        output_seqs.append(input_seqs[seq]) #  append SeqRecord to non-BWA&BLAT-aligned list.

        # write to .fastq file every 100001 reads:
        if len(output_seqs) > 100000:
            if not output_file_made:
                output_file_made = True
                with open(output_file, "w") as outfile:
                    SeqIO.write(output_seqs, outfile, "fastq")
            else:
                with open(output_file, "a") as outfile:
                    SeqIO.write(output_seqs, outfile, "fastq")
            output_seqs = []                # Reset SeqRecord list.

    # BLAT-aligned:
    else:                                               # If read is BLAT-aligned,
        contaminat_output_seqs.append(input_seqs[seq])  #  append SeqRecord to BLAT-aligned list.

        # write to .fastq file every 100001 reads:
        if len(contaminat_output_seqs) > 100000:
            if not contaminat_output_file_made:
                contaminat_output_file_made = True
                with open(contaminat_output_file, "w") as outfile:
                    SeqIO.write(contaminat_output_seqs, outfile, "fastq")
            else:
                with open(contaminat_output_file, "a") as outfile:
                    SeqIO.write(contaminat_output_seqs, outfile, "fastq")
            contaminat_output_seqs = []     # Reset SeqRecord list.

# write remaining non-BLAT-aligned reads to .fastq file:
if len(output_seqs) > 0:
    if not output_file_made:
        output_file_made = True
        with open(output_file, "w") as outfile:
            SeqIO.write(output_seqs, outfile, "fastq")
    else:
        with open(output_file, "a") as outfile:
            SeqIO.write(output_seqs, outfile, "fastq")

# write remaining BLAT-aligned reads to .fastq file:
if len(contaminat_output_seqs) > 0:
    if not contaminat_output_file_made:
        contaminat_output_file_made = True
        with open(contaminat_output_file, "w") as outfile:
            SeqIO.write(contaminat_output_seqs, outfile, "fastq")
    else:
        with open(contaminat_output_file, "a") as outfile:
            SeqIO.write(contaminat_output_seqs, outfile, "fastq")
