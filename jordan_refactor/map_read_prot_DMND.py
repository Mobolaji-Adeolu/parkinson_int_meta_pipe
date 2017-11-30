#!/usr/bin/env python

# Now with some commenting!
# CHANGES:
# - fixed sort function to return float
# - reversed the sort direction to get largest scores at top
# - multiplied align_len*3 to convert aa->nt
# - changed align_len to float(align_len) (to make sure it's a number)
# - fixed contig check to act on current read (not previous)
# - fixed checking for duplicate reads to "if query in mapped_reads:"
# - fixed WRITE OUTPUT: non-BWA&BLAT&DMD-aligned: changed "mapped_reads.add" to "break"
# - fixed WRITE OUTPUT: BWA&BLAT&DMD-aligned: changed to only write prot once

# NOTE:
# -  Sometimes reads will be matched to multiple genes. Multiples come from:
#    (1) BWA alignment itself,
#    (2) between BWA alignement and BLAT-(from contig) alignment, and
#    (3) between BWA/BLAT alignement and DMD-(from contig) alignment.
#    BLATpp doesn't allow BLAT itself to assign reads to different genes.
#    BLATpp only takes the read<->gene match with the top score. For duplicate
#    top scores, it only takes the one that shows up first in the initial sort.

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

Prot_DB= sys.argv[1]            # INPUT: AA db used for DIAMOND alignement
contig2read_file= sys.argv[2]   # INPUT: [contigID, #reads, readIDs ...]
gene2read_file= sys.argv[3]     # INPUT: [BWA&BLAT-aligned geneID, length, #reads, readIDs ...]
                                # ->OUTPUT: [BWA&BLAT&DMD-aligned gene/protID, length, #reads, readIDs ...]
gene_file= sys.argv[4]          # INPUT: BWA&BLAT-aligned geneIDs and nt seqs (.fna; fasta-format)
prot_file= sys.argv[5]          # OUTPUT: BWA&BLAT&DMD-aligned gene/protIDs and aa seqs (.faa; fasta-format)

# make dict of contigID<->readsID(s):
contig2read_map= {}
with open(contig2read_file,"r") as mapping:
    for line in mapping:
        if len(line)>5:                             # line starts with 'NODE_'
            entry= line.strip("\n").split("\t")     # break tab-separated into list
            contig2read_map[entry[0]]= entry[2:]    # key=contigID, value=list of readID(s)

# make dict of BWA&BLAT-aligned geneID<->readID(s):
gene2read_map= {}
gene_len= {}
with open(gene2read_file,"r") as mapping:
    for line in mapping:
        if len(line)>5:                             # line at least 5 characeters?
            entry= line.split("\t")
            gene2read_map[entry[0]]= entry[3:]      # key=geneID, value=list of readID(s)
            gene_len[entry[0]]= entry[1]            # key=geneID, value=gene length; to avoid recalc later

# make dict of BWA&BLAT-aligned geneID<->seq
gene_seqs= SeqIO.index(gene_file,"fasta")           # key=geneID, value=SeqRecord

# tracking DMD-assigned:
mapped_reads= set()         # tracks DMD-assigned reads
mapped_contigs= set()       # tracks DMD-assigned contigs
prot2read_map= {}           # dict of DMD-aligned protID<->readID(s)
                            #  key=protID, value=list of readID(s)

# sort function: sort by score:
# (12th field of the .blatout file)
def sortbyscore(line):
    return line[11]

# loop over remainder (after argv[5]) in sets of 3:
# (readtypes sets: contigs, merged, unmerged1, unmerged2)
for x in range((len(sys.argv)-6)/3):
    read_file= sys.argv[3*x+6]      # INPUT: non-BWA&BLAT-aligned readIDs and seqs (.fasta)
    read_seqs= SeqIO.index(read_file, os.path.splitext(read_file)[1][1:])
                                    # dict of non-BWA&BLAT-aligned read SeqRecords: key=readID
                                    #  (second argument specifies filetype, e.g., "fasta")
    DMND_tab_file= sys.argv[3*x+7]  # INPUT: DMD-aligned readIDs (.dmdout)
    output_file= sys.argv[3*x+8]    # OUTPUT: non-BWA&BLAT&DMD-aligned readIDs and seqs (.fasta)

    #####################################
    # FUNCTION:
    # add DMD-aligned reads that meet threshold
    # to the aligned gene/protID<->readID(s) dict:
    # for duplicate matches, the read is assigned to a single proteing w highest match score.

    # call as gene_map(.dmdout file, list of unmapped readIDs)
    def gene_map(tsv,unmapped):
 
        # get info from .dmdout file:
        with open(tsv,"r") as tabfile:
            Hits= []                                # List of lists containing .dmdout fields.
            for line in tabfile:                    # In the .dmdout file:
                if len(line)<2:                     # If length of line < 2,
                    continue                        #  go to next line (restart for),
                else:                               #  or else
                    Hits.append(line.split("\t"))   #  append list of tab-delimited fields to Hits list.
        
        # Sort .dmdout list by high score:
        Sorted_Hits= sorted(Hits,key=sortbyscore,reverse=True)
        
        # DMD threshold:
        identity_cutoff= 85
        length_cutoff= 0.65
        score_cutoff= 60
        
        # loop through DMD-aligned reads/contigs:
        for line in Sorted_Hits:

            # store queryID:
            query= line[0]                  # queryID= readID/contigID
            
            # process only if queryID is thus far DMD-"novel"
            # (not already recorded as DMD-matched):
            if query in contig2read_map:    # If query is a contig (searches through keys):
                if query in mapped_contigs: #  if it's a duplicate (thus w lower score than previous)
                    continue                #  skip to next entry,
                contig= True                #  or else mark as contig and continue.
            else:                           # If query isn't contig:
                if query in mapped_reads:   #  if it's duplicate (thus w lower score than previous)
                    continue                #  skip to next entry,
                contig= False               #  or else continue.
            
            # store remaining read info:
            db_match= line[1]               # proteinID
            seq_identity= float(line[2])    # sequence identity
            align_len= 3*int(line[3])       # alignment length (aa->nt)
            score= float(line[11])          # score
            
            # test thresholds:
            if seq_identity > identity_cutoff:                            # identity
                if align_len > len(read_seqs[query].seq)*length_cutoff:   # length
                    if score > score_cutoff:                              # score

                        # RECORD alignment:
                    
                        # query aligns to protein with previous alignment(s):
                        # (just look at DMD-aligned, since BWA&BLAT are genes)
                        if db_match in prot2read_map:                           # If query matches to a previously found prot:
                            if contig:                                          # If query is a contig, then
                                mapped_contigs.add(query)                       #  mark contig as assigned by DMD,
                                for read in contig2read_map[query]:             #  take all reads making up that contig
                                    if read not in prot2read_map[db_match]:     #  not already assigned to that matched prot (by DMD)
                                        if read not in mapped_reads:            #  and not already assigned by DMD to a diff prot:
                                            prot2read_map[db_match].append(read)#  and append their readIDs to aligned prot<->read dict,
                                            mapped_reads.add(read)              #  and mark them as assigned by DMD.
                            elif not contig:                                    # If query is a read, and
                                if query not in prot2read_map[db_match]:        #  it hasn't already been assigned to that matched prot:
                                    prot2read_map[db_match].append(query)       #  append its readID to aligned prot<->read dict,
                                    mapped_reads.add(query)                     #  and mark it as assigned by DMD.
                    
                        # query matches to DMD-novel prot:
                        else:                                                   # If query matches to a novel prot:
                            if contig:                                          # If query is a contig, then
                                mapped_contigs.add(query)                       #  mark contig as assigned by DMD,
                                read_count= 0
                                for read in contig2read_map[query]:             #  take all reads making up that contig
                                    if read not in mapped_reads:                #  not already assigned by DMD to a diff prot:
                                        mapped_reads.add(read)                  #  and mark them as assigned by DMD, then
                                        read_count+= 1
                                        if read_count==1:
                                            prot2read_map[db_match]= [read]     #  add its readID to aligned prot<->read dict
                                        elif read_count>1:
                                            prot2read_map[db_match].append(read)#  or append its readID to aligned prot<->read dict.
                            elif not contig:                                    # If query is a read, then
                                prot2read_map[db_match]= [query]                #  add its readID to aligned prot<->read dict
                                mapped_reads.add(query)                         #  and mark it as assigned by DMD.
                
                        continue    # If DMD-aligned query was DMD-novel and met threshold
                                    #  (and was therefore processed), skip to next query
                                    
            unmapped.add(query)     # If DMD-aligned query was DMD-novel but failed threshold
                                    #  put the queryID back in the unmapped set.

    #####################################
    
    # for DMD-aligned reads that
    # don't meet threshold:
    unmapped_reads= set()           # for set of readIDs
    unmapped_seqs= []               # for list of SeqRecords

    # process DMD-aligned reads:
    gene_map(DMND_tab_file,unmapped_reads)
    #  Novel reads that meet threshold, append to prot2read_map (aligned gene/protID<->readID(s) dict).
    #  readIDs of reads that never meet threshold are placed in unmapped_reads.

    # WRITE OUTPUT: non-BWA&BLAT&DMD-aligned readIDs
    # and seqs (.fasta)
    for read in read_seqs:                              # Take all non-BWA&BLAT-aligned reads (input to DMD)
        if read not in unmapped_reads:                  #  that weren't already in unmapped_reads;
            for prot in prot2read_map:                  #  and if the read was never aligned to any protein
                if read in prot2read_map[prot]:         #  (check against all reads in each previously-aligned protein),
                   break
            else:
                unmapped_reads.add(read)                #  then add its readID to unmapped_reads. (for/else loop)

    for read in unmapped_reads:                         # Put corresponding SeqRecords for unmapped_reads
        unmapped_seqs.append(read_seqs[read])           #  into unmapped_seqs
    with open(output_file,"w") as outfile:
        SeqIO.write(unmapped_seqs,outfile,"fasta")      #  and write it to file.


# WRITE OUTPUT: rewrite gene2read file to include DMD-aligned:
# [BWA&BLAT&DMD-aligned geneID, length, #reads, readIDs ...]
reads_count= 0
proteins= []
with open(gene2read_file,"w") as out_map:               # Delete old gene2read_file and write a new one.

    # write genes:
    for gene in gene2read_map:                          # Take each BWA&BLAT-aligned gene and
        out_map.write(gene + "\t" + gene_len[gene] + "\t" + str(len(gene2read_map[gene])))
                                                        #  write [aligned geneID, length (in nt), #reads],
        for read in gene2read_map[gene]:
            out_map.write("\t" + read.strip("\n"))  #  [readIDs ...],
        else:
            out_map.write("\n")                     #  and a new line character.

    # write proteins:
    for record in SeqIO.parse(Prot_DB,"fasta"):         # Loop through SeqRec of all prot in PROTdb:
                                                        #  (PROTdb is needed to get the aa sequence.)
        if record.id in prot2read_map:                  #  If PROTdb prot is one of the matched proteins,
            proteins.append(record)                     #  append the SeqRec to proteins list (for next file), and
            out_map.write(record.id + "\t" + str(len(record.seq)*3) + "\t" + str(len(prot2read_map[record.id])))
                                                        #  write [aligned protID, length (in nt), #reads],
            for read in prot2read_map[record.id]:
                out_map.write("\t" + read.strip("\n"))  #  [readIDs ...],
                reads_count+= 1
            else:
                out_map.write("\n")                     #  and a new line character.

# WRITE OUTPUT: BWA&BLAT&DMD-aligned gene/protIDs and aa seqs
# (.faa; fasta-format):
genes_trans= []
for gene in gene_seqs:                                  # Take each BWA&BLAT-aligned genes
    try:
        genes_trans.append(SeqRecord(seq= gene_seqs[gene].seq.translate(stop_symbol=""), id= gene_seqs[gene].id, description= gene_seqs[gene].description))
                                                        #  and translate its SeqRecord sequence to aa.
    except:
        pass
with open(prot_file,"w") as out_prot:
    SeqIO.write(genes_trans,out_prot,"fasta")           # Write aligned gene aa seqs
    SeqIO.write(proteins,out_prot,"fasta")              #  and aligned proteins aa seqs.

# print DMD stats:
print str(reads_count) + " reads were mapped with Diamond."
print "Reads mapped to " + str(len(proteins)) + " proteins."