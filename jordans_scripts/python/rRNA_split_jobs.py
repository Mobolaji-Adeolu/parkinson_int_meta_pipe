#!/usr/bin/env python

# CHANGES:
# - removed JobID1 = sys.argv[2]
# - shortened PBS submission time


import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

Input_File = sys.argv[1]
Input_Path = os.path.dirname(Input_File)
Input_FName = os.path.basename(Input_File)

#Python = "/home/j/jparkins/mobolaji/python"
Python = "/scinet/gpc/tools/Python/Python272-shared/bin/python"
Filter_rRNA = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"

PBS_Submit = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=00:15:00
#PBS -N NAME

module load gcc intel/15.0.2 openmpi java blast extras python
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Python27/Python-2.7.12/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

Preprocess_jobs = []
for split in sorted(os.listdir(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants"))):
    if split.endswith(".fastq"):
        Split_File = os.path.splitext(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", split))[0]

        COMMANDSx = [Python + " " + Filter_rRNA + " " + Split_File + ".fastq" + " " + Split_File + "_mRNA.fastq" + " " + Split_File + "_rRNA.fastq"]

        with open(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.splitext(split)[0] + "_rRNA_Filter.pbs"), "w") as PBS_script_out:
            for line in PBS_Submit.splitlines():
                if "NAME" in line:
                    line = line.replace("NAME", os.path.splitext(split)[0] + "_rRNA_Filter")
                if "COMMANDS" in line:
                    PBS_script_out.write("\n".join(COMMANDSx))
                    break
                PBS_script_out.write(line + "\n")
        JobIDx = subprocess.check_output("ssh gpc01 " + "\"" + "cd " + os.path.dirname(Split_File) + ";" + "qsub" + " " + os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.splitext(split)[0] + "_rRNA_Filter.pbs") + "\"", shell=True)
        Preprocess_jobs.append(JobIDx.strip("\n"))

for split in sorted(os.listdir(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants"))):
    if split.split("_paired_n_contaminants_split_")[0].endswith("1"):
        Split_File1 = os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", split.split("_paired_n_contaminants_split_")[0][:-1] + "1" + "_paired_n_contaminants_split_" + split.split("_paired_n_contaminants_split_")[1])
        Split_File2 = os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", split.split("_paired_n_contaminants_split_")[0][:-1] + "2" + "_paired_n_contaminants_split_" + split.split("_paired_n_contaminants_split_")[1])

        COMMANDSy = [Python + " " + Filter_rRNA + " " + Split_File1 + " " + os.path.splitext(Split_File1)[0] + "_mRNA.fastq" + " " + os.path.splitext(Split_File1)[0] + "_rRNA.fastq" + " " + Split_File2 + " " + os.path.splitext(Split_File2)[0] + "_mRNA.fastq" + " " + os.path.splitext(Split_File2)[0] + "_rRNA.fastq"]

        with open(os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.splitext(split)[0] + "_rRNA_Filter.pbs"), "w") as PBS_script_out:
            for line in PBS_Submit.splitlines():
                if "NAME" in line:
                    line = line.replace("NAME", os.path.splitext(split)[0] + "_rRNA_Filter")
                if "COMMANDS" in line:
                    PBS_script_out.write("\n".join(COMMANDSy))
                    break
                PBS_script_out.write(line + "\n")
        JobIDy = subprocess.check_output("ssh gpc01 " + "\"" + "cd " + os.path.dirname(Split_File1) + ";" + "qsub" + " " + os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.splitext(split)[0] + "_rRNA_Filter.pbs") + "\"", shell=True)
        Preprocess_jobs.append(JobIDy.strip("\n"))

print ":".join(Preprocess_jobs[-10:])
