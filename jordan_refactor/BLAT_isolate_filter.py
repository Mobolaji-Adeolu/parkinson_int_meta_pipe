#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - closed the blatout_tab_file, after blatout_readID is made

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

# INPUTS:
BWA_file = sys.argv[1]          # INPUT: BWA-aligned reads (.fastq)
BWA_seqs = SeqIO.to_dict(SeqIO.parse(BWA_file, "fastq"))    # dict of BWA-aligned read SeqRecords
input_file = sys.argv[2]        # INPUT: non-BWA-aligned reads (.fastq); input to BLAT
input_seqs = SeqIO.to_dict(SeqIO.parse(input_file, "fastq"))# dict of non-BWA-aligned read SeqRecords
                                                            #  SeqIO.parse parses the .fastq one read at a time, 
                                                            #  creates a SeqRecord object > SeqIO.to_dict places 
                                                            #  SeqRecord into dictionary using readID as dict key.
blatout_tab_file = sys.argv[3]  # INPUT: BLAT output (.blatout) (BLAT uses non-BWA aligned reads as input)
                                #  ****IMPORTANT**** .blatout has identical readIDs on ajacent lines.

# OUTPUTS:
nonaligned_file = sys.argv[4]   # OUTPUT: non-BWA&BLAT-aligned reads (.fastq)
nonaligned_seqs = []            # list of non-BWA&BLAT-aligned read SeqRecords
nonaligned_file_made = False
BLAT_file = sys.argv[5]         # OUTPUT: BLAT-aligned reads (.fastq)
BLAT_seqs = []                  # list of BLAT-aligned read SeqRecords
BLAT_file_made = False
BWABLAT_file = sys.argv[6]      # OUTPUT: BWA&BLAT-aligned reads (.fastq)
BWABLAT_seqs = []               # liast of BWA&BLAT-aligned read SeqRecords

# copy all BWA-aligned reads into 
# BWA&BLAT-aligned reads output file:
shutil.copyfile(BWA_file, BWABLAT_file)

# make a list of BLAT-aligned readIDs:
with open(blatout_tab_file, "r") as tabfile:    # Open .blatout file.
    blatout_readID = []
    query_seq = ""
    for line in tabfile:                        # In the .blatout file:
        if len(line) < 2:                       # If length of line < 2,
            continue                            #  go to next line (restart for).
        line_parts = line.split("\t")           # Split line by tab separator (Element 1=readID).
        if query_seq == line_parts[0]:          # If readID is empty or same as last line,
            continue                            #  skip to next line.
        else:
            query_seq = line_parts[0]           # If readID exists, append to a
            blatout_readID.append(query_seq)    #  list of BLAT-aligned readIDs. 

# separate non-BWA-aliged reads
# according to BLAT-alignment:
for readID in input_seqs:                           # Look though all non-BWA-aligned readIDs:

    # non-BLAT-aligned:                      
    if readID not in blatout_readID:                # If read is non-BLAT-aligned,  
        nonaligned_seqs.append(input_seqs[readID])  #  append SeqRecord to non-BWA&BLAT-aligned list.

        # write to .fastq file every 100001 reads:
        if len(nonaligned_seqs) > 100000:
            if not nonaligned_file_made:
                nonaligned_file_made = True
                with open(nonaligned_file, "w") as outfile:
                    SeqIO.write(nonaligned_seqs, outfile, "fastq")
            else:
                with open(nonaligned_file, "a") as outfile:
                    SeqIO.write(nonaligned_seqs, outfile, "fastq")
            nonaligned_seqs = []                    # Reset SeqRecord list.

    # BLAT-aligned:
    else:                                           # If read is BLAT-aligned,
        BLAT_seqs.append(input_seqs[readID])        #  append SeqRecord to BLAT-aligned list,
        BWABLAT_seqs.append(input_seqs[readID])     #  and to the BWA&BLAT-aligned list.

        # write to .fastq files every 100001 reads:
        if len(BLAT_seqs) > 100000:
            if not BLAT_file_made:
                BLAT_file_made = True
                with open(BLAT_file, "w") as outfile:
                    SeqIO.write(BLAT_seqs, outfile, "fastq")
            else:
                with open(BLAT_file, "a") as outfile:
                    SeqIO.write(BLAT_seqs, outfile, "fastq")

            with open(BWABLAT_file, "a") as outfile:
                    SeqIO.write(BWABLAT_seqs, outfile, "fastq")
            BLAT_seqs = []                          # Reset SeqRecord list.
            BWABLAT_seqs = []                       # Reset SeqRecord list.

# write remaining non-BLAT-aligned reads to .fastq file:
if len(nonaligned_seqs) > 0:
    if not nonaligned_file_made:
        nonaligned_file_made = True
        with open(nonaligned_file, "w") as outfile:
            SeqIO.write(nonaligned_seqs, outfile, "fastq")
    else:
        with open(nonaligned_file, "a") as outfile:
            SeqIO.write(nonaligned_seqs, outfile, "fastq")

# write remaining BLAT-aligned reads to .fastq files:
if len(BLAT_seqs) > 0:
    if not BLAT_file_made:
        BLAT_file_made = True
        with open(BLAT_file, "w") as outfile:
            SeqIO.write(BLAT_seqs, outfile, "fastq")
    else:
        with open(BLAT_file, "a") as outfile:
            SeqIO.write(BLAT_seqs, outfile, "fastq")

    with open(BWABLAT_file, "a") as outfile:
            SeqIO.write(BWABLAT_seqs, outfile, "fastq")
