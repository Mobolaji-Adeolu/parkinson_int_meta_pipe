#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
import multiprocessing
from Bio import SeqIO

splitlen = int(sys.argv[1])
input_file = sys.argv[2]
output_folder = sys.argv[3]

def check_file(ending, count):
    output_file = os.path.join(output_folder, os.path.splitext(os.path.basename(input_file))[0] + "_split_" + str(count).zfill(3) + ending)
    alt_output_file = os.path.join(output_folder, os.path.splitext(os.path.basename(input_file))[0][:-1] + "_split_" + str(count).zfill(3) + "_" + os.path.splitext(os.path.basename(input_file))[0][-1:] + ending)
    result = False
    if os.path.exists(output_file) or os.path.exists(alt_output_file):
        result = True
    return result

def write_fastq(fastq_dict, count, mode):
    output_file = os.path.join(output_folder, os.path.splitext(os.path.basename(input_file))[0] + "_split_" + str(count).zfill(3) + ".fastq")
    with open(output_file, mode) as outfile:
        for seq_id in fastq_dict:
            outfile.write("@" + seq_id + "\n")
            outfile.write(fastq_dict[seq_id][0] + "\n")
            outfile.write("+" + "\n")
            outfile.write(fastq_dict[seq_id][1] + "\n")

def write_fasta(fasta_dict, count, mode):
    output_file = os.path.join(output_folder, os.path.splitext(os.path.basename(input_file))[0] + "_split_" + str(count).zfill(3) + ".fasta")
    with open(output_file, mode) as outfile:
        for seq_id in fasta_dict:
            outfile.write(">" + seq_id + "\n")
            outfile.write(fasta_dict[seq_id][0] + "\n")

os.chdir(os.path.dirname(input_file))

if input_file.endswith("q"):
    read_count = 0
    with open(input_file, "r") as infile:
        for line in infile:
            read_count += 1
        read_count = read_count/4
    if read_count % splitlen == 0:
        File_count = read_count / splitlen # added such that file with 0 reads isn't created
    else:
        File_count = read_count / splitlen + 1
    for n in range(File_count):
        Start = n*splitlen
        m = n + 1
        Stop = m*splitlen
        Temp_sequences = {}
        if check_file(".fastq", m + 1):
            continue
        with open(input_file, "r") as infile:
            seq_id = ""
            seq_count = 0
            line_count = 0
            for line in infile:
                line_class = line_count % 4
                if seq_count < Start:
                    if line_class == 0:
                        seq_count += 1
                elif seq_count >= Start and seq_count < Stop:
                    if line_class == 0 and line.startswith("@"):
                        seq_id = line[1:].split(" ")[0].strip("\n")
                        Temp_sequences[seq_id] = ["",""]
                    elif seq_id != "":
                        if line_class == 1 and Temp_sequences[seq_id][0] == "":
                            Temp_sequences[seq_id][0] = line.strip("\n")
                        elif line_class == 3 and Temp_sequences[seq_id][1] == "":
                            Temp_sequences[seq_id][1] = line.strip("\n")
                            seq_count += 1
                elif seq_count == Stop:
                    break
                line_count += 1
        write_fastq(Temp_sequences, m, "w")

elif input_file.endswith("a"):
    prot_count = 0
    with open(input_file, "r") as infile:
        for line in infile:
            if line.startswith(">"):
                prot_count += 1
    if prot_count % splitlen == 0:
        File_count = prot_count / splitlen # added such that file with 0 reads isn't created
    else:
        File_count = prot_count / splitlen + 1
    for n in range(File_count):
        Start = n*splitlen
        m = n + 1
        Stop = m*splitlen
        Temp_sequences = {}
        if check_file(".fasta", m + 1):
            continue
        with open(input_file, "r") as infile:
            seq_id = ""
            seq_count = 0
            for line in infile:
                if line.startswith(">"):
                    if seq_count < Start:
                        seq_count += 1
                    elif seq_count >= Start and seq_count < Stop:
                        seq_id = line[1:].split(" ")[0].strip("\n")
                        Temp_sequences[seq_id] = [""]
                        seq_count += 1
                    elif seq_count == Stop:
                        break
                elif seq_id != "":
                    Temp_sequences[seq_id][0] += line.strip("\n")
        write_fasta(Temp_sequences, m, "w")
else:
    print "Unrecognized File Type"
