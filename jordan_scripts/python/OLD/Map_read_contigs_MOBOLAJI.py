#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

paired_file1 = sys.argv[1]
paired1_seqs = SeqIO.index(paired_file1, "fastq")
paired_file2 = sys.argv[2]
paired2_seqs = SeqIO.index(paired_file2, "fastq")
unpaired_file = sys.argv[3]
unpaired_seqs = SeqIO.index(unpaired_file, "fastq")
paired_sam_file = sys.argv[4]
unpaired_sam_file = sys.argv[5]
output_file = sys.argv[6]


contig2read_map = {}
unmapped_paired_reads = set()
unmapped_paired1_seqs = []
unmapped_paired2_seqs = []
unmapped_unpaired_reads = set()
unmapped_unpaired_seqs = []


def contig_map(sam, unmapped):
    with open(sam, "r") as samfile:
        for line in samfile:
            if not line.startswith("@") and len(line) > 1:
                line_parts = line.split("\t")
                flag = bin(int(line_parts[1]))[2:].zfill(11)
                if flag[8] == "0":
                    if line_parts[2] in contig2read_map:
                        if line_parts[0] not in contig2read_map[line_parts[2]]:
                            contig2read_map[line_parts[2]].append(line_parts[0])
                    else:
                        contig2read_map[line_parts[2]] = [line_parts[0]]
                elif flag[8] == "1":
                    unmapped.add(line_parts[0])

contig_map(paired_sam_file, unmapped_paired_reads)
contig_map(unpaired_sam_file, unmapped_unpaired_reads)

for read in unmapped_paired_reads:
    unmapped_paired1_seqs.append(paired1_seqs[read])
    unmapped_paired2_seqs.append(paired2_seqs[read])
with open(os.path.splitext(paired_file1)[0] + "_unmapped.fastq", "w") as out:
    SeqIO.write(unmapped_paired1_seqs, out, "fastq")
with open(os.path.splitext(paired_file2)[0] + "_unmapped.fastq", "w") as out:
    SeqIO.write(unmapped_paired2_seqs, out, "fastq")

for read in unmapped_unpaired_reads:
    unmapped_unpaired_seqs.append(unpaired_seqs[read])
with open(os.path.splitext(unpaired_file)[0] + "_unmapped.fastq", "w") as out:
    SeqIO.write(unmapped_unpaired_seqs, out, "fastq")

with open(output_file, "w") as outfile:
    for contig in contig2read_map:
        outfile.write(contig + "\t" + str(len(contig2read_map[contig])))
        for read in contig2read_map[contig]:
            outfile.write("\t" + read)
        else:
            outfile.write("\n")
