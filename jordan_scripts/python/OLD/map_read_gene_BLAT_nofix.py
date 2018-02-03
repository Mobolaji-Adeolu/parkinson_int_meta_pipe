#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - changed align_len to float(align_len) (to make sure it's a number)

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

DNA_DB = sys.argv[1]                # INPUT: DNA db used for BLAT alignement
contig2read_file = sys.argv[2]      # INPUT: [contigID, #reads, readIDs ...]
gene2read_file = sys.argv[3]        # INPUT: [BWA-aligned geneID, length, #reads, readIDs ...]
                                    # OUTPUT: [BWA&BLAT-aligned geneID, length, #reads, readIDs ...]
gene_file = sys.argv[4]             # OUTPUT: BWA&BLAT-aligned geneIDs and seqs (.fna; fasta-format)

# make dict of contigID<->readsID(s):
contig2read_map = {}
with open(contig2read_file, "r") as mapping:
    for line in mapping:
        if len(line) > 5:                           # line starts with 'NODE_'
            entry = line.strip("\n").split("\t")    # break tab-separated into list
            contig2read_map[entry[0]] = entry[2:]   # key=contigID, value=list of readID(s)

# make dict of BWA-aligned geneID<->readID(s):
gene2read_map = {}
with open(gene2read_file, "r") as mapping:
    for line in mapping:
        if len(line) > 5:                           # line at least 5 characeters?
            entry = line.split("\t")
            gene2read_map[entry[0]] = entry[3:]     # key=geneID, value=list of readID(s)


mapped_reads = set()                # tracks BLAT-assigned reads

# sort function: sort by score:
# (12th field of the .blatout file)
def sortbyscore(line):
    return line[11]

# loop over remainder (after argv[4]) in sets of 3:
# (readtypes sets: contigs, merged, unmerged1, unmerged2)
for x in range((len(sys.argv) - 5) / 3):
    read_file = sys.argv[3 * x + 5]     # INPUT: non-BWA-aligned readIDs and seqs (.fasta)
    read_seqs = SeqIO.index(read_file, os.path.splitext(read_file)[1][1:])
                                        # dict of non-BWA-aligned read SeqRecords: key=readID
                                        #  (second argument specifies filetype, e.g., "fasta")
    BLAT_tab_file = sys.argv[3 * x + 6] # INPUT: BLAT-aligned readIDs (.blatout)
                                        #  Note: identical readIDs on ajacent lines.
    output_file = sys.argv[3 * x + 7]   # OUTPUT: non-BWA&BLAT-aligned readIDs and seqs (.fasta)

    #####################################
    # FUNCTION:
    # add BLAT-aligned reads that meet threshold
    # to the aligned geneID<->readID(s) dict:
    # for duplicate matches, the read is assigned to a single gene with highest match score.
    
    # for gene_map(.blatout file, list of unmapped readID)
    def gene_map(tsv, unmapped):
    
        # get info from .blatout file:
        with open(tsv, "r") as tabfile:
            Hits = []                                       # List of .blatout fields.
            for line in tabfile:                            # In the .blatout file:
                if len(line) < 2:                           # If length of line < 2,
                    continue                                #  go to next line (restart for),
                else:                                       #  or else
                    Hits.append(line.split("\t"))           #  append list of tab-delimited fields to Hits list.
    
        Sorted_Hits = sorted(Hits, key = sortbyscore)       # Sort .blatout list by score.

        # BLAT threshold:
        identity_cutoff = 85
        length_cutoff = 0.65
        score_cutoff = 60

        # loop through BLAT-aligned reads:
        query = ""                          # Initial "previous" readID.
        for line in Sorted_Hits:
            if query in contig2read_map:    # Check if previous read is a contig itself (searches keys).
                contig = True
            else:
                contig = False
            if query == line[0]:            # If readID is empty or same as previous readID
                continue                    #  (which would have a higher score), skip to next entry.
            else:
                query = line[0]             # readID
                db_match = line[1]          # geneID
                seq_identity = line[2]      # sequence identity
                align_len = line[3]         # alignment length
                score = line[11]            # score
            
            # test thresholds:
            if float(seq_identity) > float(identity_cutoff):                        #  identity
                if float(align_len) > len(read_seqs[query].seq) * length_cutoff:    #  length
                    if float(score) > float(score_cutoff):                          #  score
                    
                        # read matches to geneID already found:
                        if db_match in gene2read_map:                           # If read matches to a previously found gene:
                            if contig:                                          # If read is a contig, then
                                for read in contig2read_map[query]:             #  take all reads making up that contig
                                    if read not in gene2read_map[db_match]:     #  that aren't already assigned to that matched gene,
                                        if read not in mapped_reads:            #  and haven't already been assigned by BLAT to diff gene:
                                            gene2read_map[db_match].append(read)#  and append their readIDs to the aligned gene<->read dict,
                                            mapped_reads.add(read)              #  and mark them as assigned by BLAT.
                        
                                                                                # **** Checking to see if previousy found by BLAT is only
                                                                                #  for reads in a contig.........
                        
                            elif not contig:                                    # If read isn't a contig, and
                                if query not in gene2read_map[db_match]:        #  it hasn't already been assigned to that matched gene:
                                    gene2read_map[db_match].append(query)       #  append its readID to the aligned gene<->read dict,
                                    mapped_reads.add(query)                     #  and mark it as assigned by BLAT.
                    
                        # read matches to novel geneID:
                        else:                                                   # If read matches to a novel gene:
                            if contig:                                          # If read is a contig, then
                                read_count = 0
                                for read in contig2read_map[query]:             #  take all reads making up that contig
                                    if read not in mapped_reads:                #  and haven't already been assigned by BLAT to diff gene:
                                        mapped_reads.add(read)                  #  and mark them as assigned by BLAT, then
                                        read_count += 1
                                        if read_count == 1:
                                            gene2read_map[db_match] = [read]    #  add its readID to the aligned gene<->read dict
                                        elif read_count > 1:
                                            gene2read_map[db_match].append(read)#  or append its readID to the aligned gene<->read dict.
                            elif not contig:                                    # If read isn't a contig, then
                                gene2read_map[db_match] = [query]               #  add its readID to the aligned gene<->read dict
                                mapped_reads.add(query)                         #  and mark it as assigned by BLAT.
                                
                        continue    # If the BLAT-aligned read meets threshold and is kept,
                                    #  skip to the next BLAT-aligned read.
                                    
                                    
            unmapped.add(query)     # If BLAT-aligned read was NOT same as a previous entry AND
                                    #  failed the threshold, put the readID in the unmapped set.

    #####################################

    # for BLAT-aligned reads that don't meet threshold:
    unmapped_reads = set()          # for set of readIDs
    unmapped_seqs = []              # for list of SeqRecords

    # process BLAT-aligned reads:
    gene_map(BLAT_tab_file, unmapped_reads)
    #  unique reads that meet threshold, append to gene2read_map (aligned geneID<->readID(s) dict)
    #  readIDs of reads that never meet threshold are placed in unmapped_reads


    # WRITE OUTPUT: non-BWA&BLAT-aligned readIDs and seqs (.fasta)
    for read in read_seqs:                              # Take all non-BWA-aligned reads (input to BLAT)
        if read not in unmapped_reads:                  #  that weren't already in unmapped_reads;
        
            # for/else loop:
            for gene in gene2read_map:                  #  and if the read was not aligned to any gene
                if read in gene2read_map[gene]:         #  (check against all reads in each previously-aligned gene),
                    break
            else:
                unmapped_reads.add(read)                #  then add it's readID to unmapped_reads.

    for read in unmapped_reads:                         # Put corresponding SeqRecords for unmapped_reads
        unmapped_seqs.append(read_seqs[read])           #  into unmapped_seqs

    with open(output_file, "w") as outfile:
        SeqIO.write(unmapped_seqs, outfile, "fasta")    #  and write it to file.



# REWRITE gene2read OUTPUT file to include BLAT-aligned:
# [BWA&BLAT-aligned geneID, length, #reads, readIDs ...]
reads_count = 0
genes = []
with open(gene2read_file, "w") as out_map:              # Delete old gene2read_file and writes a new one.
    for record in SeqIO.parse(DNA_DB, "fasta"):         # Loop through SeqRec of all genes in DNA db:
                                                        #  (DNA db is needed to get the sequence.)
        if record.id in gene2read_map:                  # If DNA db gene is one of the matched genes,
            genes.append(record)                        #  append the SeqRec to genes.
            out_map.write(record.id + "\t" + str(len(record.seq)) + "\t" + str(len(gene2read_map[record.id])))
                                                        # Write [aligned geneID, length, #reads],
            for read in gene2read_map[record.id]:
                out_map.write("\t" + read.strip("\n"))  #  then write [readIDs ...]
                reads_count += 1
            out_map.write("\n")                         #  and a new line character.

# WRITE OUTPUT: BWA&BLAT-aligned geneIDs and seqs (.fna; fasta-format)
with open(gene_file, "w") as outfile:
    SeqIO.write(genes, outfile, "fasta")

print str(reads_count) + " reads were mapped with BWA and BLAT"
print "Reads mapped to " + str(len(genes)) + " genes."