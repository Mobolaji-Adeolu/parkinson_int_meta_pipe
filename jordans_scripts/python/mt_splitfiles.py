#!/usr/bin/env python

# CHANGED:
# - now creates "splitfiles" even if num reads <= split number exactly

import sys
import os
import os.path
import subprocess
import multiprocessing

#Python = "/home/j/jparkins/mobolaji/python"
Python = "/scinet/gpc/tools/Python/Python272-shared/bin/python"
#File_splitter = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/File_splitter.py"
File_splitter = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/file_splitter.py"
# Sort_Reads = "/home/j/jparkins/mobolaji/Read_Classification/Sort_Reads.py"
Sort_Reads = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/mt_sortreads.py"

input_folder = sys.argv[1]
splitreads = int(sys.argv[2])
splitreads_min = 100000

print "Remember to run " + Sort_Reads + " on your reads before splitting."

for genome in sorted(os.listdir(input_folder)):
    print genome
    if genome.endswith("1.fastq"):
        genome_name = os.path.splitext(genome)[0][:-1]
        line_count = 0
        with open(os.path.join(input_folder, genome), "r") as counting_file:
            for line in counting_file:
                line_count += 1
        if line_count/4 >= splitreads or line_count/4 < splitreads: # always make splitfiles now
                                                                    # even if num reads < split number
            try:
                os.mkdir(os.path.join(input_folder, genome_name))
            except:
                pass
            subprocess.call([Python, File_splitter, str(splitreads), str(splitreads_min), os.path.join(input_folder, genome), os.path.join(input_folder, genome_name)])
            subprocess.call([Python, File_splitter, str(splitreads), str(splitreads_min), os.path.join(input_folder, genome_name + "2.fastq"), os.path.join(input_folder, genome_name)])

            # move '_1' and '_2' to just before '.fastq':
            for genome_split in os.listdir(os.path.join(input_folder, genome_name)):
                if genome_split.endswith("_1.fastq"):
                    name = os.path.join(input_folder,  genome_name, genome_split[:-7])
                    if name not in file_list:
                        file_list.append(name)
                elif genome_split.endswith("_2.fastq"):
                    continue
                else:
                    try:
                        prefix = genome_split.split("_split_")[0]
                        renamed = prefix[:-1] + "_split_" + genome_split.split("_split_")[1]
                        renamed_full = os.path.join(input_folder,  genome_name, os.path.splitext(renamed)[0] + "_")
                        if prefix.endswith("1"):
                            os.rename(os.path.join(input_folder, genome_name, genome_split), renamed_full + "1.fastq")
                            if renamed_full not in file_list:
                                file_list.append(os.path.join(input_folder,  genome_name, os.path.splitext(renamed)[0] + "_"))
                        elif prefix.endswith("2"):
                            os.rename(os.path.join(input_folder, genome_name, genome_split), renamed_full + "2.fastq")
                    except:
                        pass
