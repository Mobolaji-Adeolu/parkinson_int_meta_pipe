#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
import multiprocessing
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

input_file= sys.argv[1]
output_file= sys.argv[2]

seqs= {}

with open(input_file,"r") as infile:
    seq_found = False
    seq_id = ""
    skip_plus = True
    for line in infile:
        if line.startswith("@") and not seq_found:
        
            #####################################
            # EXTRACTING SEQ_ID
#            seq_id = line[1:].split(" ")[0].strip("\n")            # choosing 1st string 1st space
            seq_id = ":".join(line.strip("\n").split(":")[1:-2])    # '@101115_D35P56SR:1:1:11664:1022:0:1'
                                                                    # -> '1:1:11664:1022'
            #####################################
            
            seqs[seq_id] = ["",""]
            seq_found = True
        elif seq_found and seqs[seq_id][0] == "":
            seqs[seq_id][0] = line.strip("\n")
        elif line.startswith("+") and skip_plus:
            skip_plus = False
        elif seq_found and not skip_plus and seqs[seq_id][0] != "":
            seqs[seq_id][1] = line.strip("\n")
            seq_found = False
            skip_plus = True

def write_fastq(fastq_dict, mode):
    with open(output_file, mode) as outfile:
        for seq_id in fastq_dict:
            outfile.write("@" + seq_id + "\n")
            outfile.write(fastq_dict[seq_id][0] + "\n")
            outfile.write("+" + "\n")
            outfile.write(fastq_dict[seq_id][1] + "\n")

file_made = False
sorted_seqs = sorted(seqs)
seq_records = {}

for seq_id in sorted_seqs:
    if seq_id.endswith(":1") or seq_id.endswith(":2") or seq_id.endswith(":3") or seq_id.endswith(":4"):
        seq_records[seq_id[:-2]] = seqs[seq_id]
    else:
        seq_records[seq_id] = seqs[seq_id]
    if len(seq_records) == 100000:
        if not file_made:
            file_made = True
            write_fastq(seq_records, "w")
        else:
            write_fastq(seq_records, "a")
        seq_records = {}
if len(seq_records) > 0:
    if not file_made:
        file_made = True
        write_fastq(seq_records, "w")
    else:
        write_fastq(seq_records, "a")