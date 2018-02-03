#!/usr/bin/env python

# CHANGES:
# - removed JobID = sys.argv[4]
# - shortened PBS submission time

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

Input_Path = sys.argv[1]
Output_Path = sys.argv[2]
Threads = sys.argv[3]

#Python = "/home/j/jparkins/mobolaji/python"
Python = "/scinet/gpc/tools/Python/Python272-shared/bin/python"
Detect = "/home/j/jparkins/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"

PBS_Submit = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=12:00:00
#PBS -N NAME

module load gcc intel/15.0.2 openmpi java blast extras python
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Python27/Python-2.7.12/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

Detect_jobs = []
for split in sorted(os.listdir(Input_Path)):
    if split.endswith(".fasta"):
        Split_File = os.path.join(Input_Path, split)
        Split_Name = os.path.splitext(split)[0]
        Split_Folder = os.path.join(Output_Path, "Detect")
        Split_Out = os.path.join(Split_Folder, Split_Name + ".detect")
        Split_Top = os.path.join(Split_Folder, Split_Name + ".toppred")

        COMMANDSx = ["mkdir -p " + Split_Folder,
        Python + " " + Detect + " " + Split_File + " " + "--output_file" + " " + Split_Out + " " + "--top_predictions_file" + " " + Split_Top + " " + "--num_threads" + " " + Threads
        ]

        with open(os.path.join(Input_Path, Split_Name + "_Detect_split.pbs"), "w") as PBS_script_out:
            for line in PBS_Submit.splitlines():
                if "NAME" in line:
                    line = line.replace("NAME", Split_Name + "_Detect_split")
                if "COMMANDS" in line:
                    PBS_script_out.write("\n".join(COMMANDSx))
                    break
                PBS_script_out.write(line + "\n")
        JobIDx = subprocess.check_output("ssh gpc01 " + "\"" + "cd " + Input_Path + ";" + "qsub" + " " + os.path.join(Input_Path, Split_Name + "_Detect_split.pbs") + "\"", shell=True)
        Detect_jobs.append(JobIDx.strip("\n"))

print ":".join(Detect_jobs[-10:])