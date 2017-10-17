#!/usr/bin/env python

# Oct 17, 2017
# ---------------------------------------
# This appears to be the main file for the entire pipeline

import sys
import os
import os.path
import subprocess
import multiprocessing

# this part's too messy.  we can clean this up
cdhit_dup = "/home/j/jparkins/mobolaji/Tools/CDHIT/cd-hit-v4.6.6-2016-0711/cd-hit-auxtools/cd-hit-dup"
Timmomatic = "/home/j/jparkins/mobolaji/Tools/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar"
AdapterRemoval = "/home/j/jparkins/mobolaji/Tools/AdapterRemoval/adapterremoval-2.2.0/build/AdapterRemoval"
vsearch = "/home/j/jparkins/mobolaji/Tools/vsearch/vsearch-2.4.2-linux-x86_64/bin/vsearch"
UniVec_Core = "/home/j/jparkins/mobolaji/Databases/UniVec_Core.fasta"
Flash = "/home/j/jparkins/mobolaji/Tools/Flash/FLASH-1.2.11/flash"
Perl = "/home/j/jparkins/mobolaji/perl"
Perl_Script_Dir = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Xuejian"
Python = "/home/j/jparkins/mobolaji/python"
BWA = "/home/j/jparkins/mobolaji/Tools/BWA/bwa-0.7.5a/bwa"
SAMTOOLS = "/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/samtools"
BLAT = "/home/j/jparkins/mobolaji/Tools/pBLAT/pblat/pblat"
DIAMOND = "/home/j/jparkins/mobolaji/Tools/Diamond/diamond"
Blastdbcmd = "/home/j/jparkins/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/blastdbcmd"
Makeblastdb = "/home/j/jparkins/mobolaji/Tools/BLAST+/ncbi-blast-2.5.0+/bin/makeblastdb"
DNA_DB = "/scratch/j/jparkins/mobolaji/Microbial_cds_db/microbial_all_cds.fasta"
DNA_DB_Prefix = os.path.splitext(DNA_DB)[0]
DNA_DB_Extension = os.path.splitext(DNA_DB)[1]
Prot_DB = "/scratch/j/jparkins/mobolaji/NCBI_nr_db/nr"
Host = "/home/j/jparkins/mobolaji/Databases/Mouse_cds.fasta"
BBMap_Dir = "/home/j/jparkins/mobolaji/Tools/BBMap/bbmap"
Fastqc = "/home/j/jparkins/mobolaji/Tools/FastQC/fastqc"
Adapter = "/home/j/jparkins/mobolaji/Tools/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa"
Infernal = "/home/j/jparkins/mobolaji/Tools/Infernal/infernal-1.1.2-linux-intel-gcc/binaries/cmscan"
Rfam = "/home/j/jparkins/mobolaji/Databases/Rfam_rRNA.cm"
Filter_rRNA = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Filter.py"
Reduplicate = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Reduplicate.py"
Map_reads_contigs = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Map_read_contigs.py"
Paired_Reads_Filter = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Paired_Reads_Filter.py"
BLAT_Contaminant_Filter = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/BLAT_Contaminant_Filter.py"
File_splitter = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/File_splitter.py"
Sort_Reads = "/home/j/jparkins/mobolaji/Read_Classification/Sort_Reads.py"
rRNA_Split_Jobs = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Split_Jobs.py"
Map_reads_gene_BWA = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Map_read_gene_BWA.py"
Map_reads_gene_BLAT = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Map_read_gene_BLAT.py"
Map_reads_prot_DMND = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Map_read_prot_DMND.py"
Spades = "/home/j/jparkins/mobolaji/Tools/SPAdes/SPAdes-3.9.1-Linux/bin/spades.py"

EC_Annotation_Prep = "/home/j/jparkins/mobolaji/EC_Prediction_Scripts/0_Preprocess_Input.py"
Detect_Submit = "/home/j/jparkins/mobolaji/EC_Prediction_Scripts/1-1a_Detect_Submission.py"
EC_Annotation_Post = "/home/j/jparkins/mobolaji/EC_Prediction_Scripts/4a_EC_Consolidation.py"
Detect = "/home/j/jparkins/mobolaji/Tools/UpdatedDETECT_V2.0/detect_leon.py"
Priam = "/home/j/jparkins/mobolaji/Tools/PRIAM/PRIAM_search.jar"
BLAST_dir = "/home/j/jparkins/mobolaji/Tools/BLAST/blast-2.2.26/bin/"
SWISS_PROT = "/home/j/jparkins/mobolaji/Databases/uniprot_sprot_annotated.fasta"

Nodes = "/home/j/jparkins/mobolaji/Databases/taxdump/nodes.dmp"
Names = "/home/j/jparkins/mobolaji/Databases/taxdump/names.dmp"
Annotated_taxid = "/home/j/jparkins/mobolaji/Read_Classification/Get_TaxID.py"
Contrain_classification = "/home/j/jparkins/mobolaji/Read_Classification/Constrain_Classification.py"
Classification_combine = "/home/j/jparkins/mobolaji/Read_Classification/Combine_WEVOTE.py"
WEVOTE = "/home/j/jparkins/mobolaji/Tools/WEVOTE/WEVOTE/run_WEVOTE_PIPELINE.sh"
WEVOTEDB = "/home/j/jparkins/mobolaji/Tools/WEVOTE/WEVOTE/WEVOTEDB"
accession2taxid = "/scratch/j/jparkins/mobolaji/accession2taxid/accession2taxid"
Kaiju = "/home/j/jparkins/mobolaji/Tools/Kaiju/kaiju-v1.4.5-linux-x86_64-static/bin/kaiju"
Kaiju2krona = "/home/j/jparkins/mobolaji/Tools/Kaiju/kaiju-v1.4.5-linux-x86_64-static/bin/kaiju2krona"
ktImportText = "/home/j/jparkins/mobolaji/Tools/Krona/KronaTools-2.7/bin/ktImportText"
Centrifuge = "/home/j/jparkins/mobolaji/Tools/Centrifuge/centrifuge/centrifuge"
Centrifuge_report = "/home/j/jparkins/mobolaji/Tools/Centrifuge/centrifuge/centrifuge-kreport"
kSLAM = "/home/j/jparkins/mobolaji/Tools/k-SLAM/k-SLAM/SLAM"

RPKM = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/RPKM.py"

Threads = str(multiprocessing.cpu_count())

# These are all different types of the same call to SciNet, using different configs
# clean it up

PBS_Submit_LowMem = """#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=12:00:00
#PBS -N NAME

module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras python
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Python27/Python-2.7.12/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

PBS_Submit_HighMem = """#!/bin/bash
#PBS -l nodes=1:m64g:ppn=16,walltime=12:00:00 -q sandy
#PBS -N NAME

module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras python
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=16
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Python27/Python-2.7.12/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

PBS_Submit_vHighMem = """#!/bin/bash
#PBS -l nodes=1:m128g:ppn=16,walltime=1:00:00 -q sandy
#PBS -N NAME

module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras python
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=20
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Python27/Python-2.7.12/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""

input_folder = sys.argv[1]
output_folder = sys.argv[2]

print "Remember to run " + Sort_Reads + " on your reads before running this pipeline."

# As the pipeline says, this is not the first stage.
# The first stage is "Sort_Reads"
# and something to generate fastq, and also splits them

for genome in sorted(os.listdir(input_folder)):
	# one main for-loop.  Script seems to loop through all the files in a dir.
    file_list = []
    Network_list = []
    if genome.endswith("1.fastq"):
		# this part splits the fastq file into segments.  why in such an odd way?
		# also recall that the fastq format has artifacts we need to discard 
		# this step could be called the "fastq cleanup"
		# recall that there was also a claim that the program fails on sufficiently small datasets... (WTF??)
		# I guess the important part was getting the main procedures outlined.
        genome_name = os.path.splitext(genome)[0][:-1]
        line_count = 0
        with open(os.path.join(input_folder, genome), "r") as counting_file:
            for line in counting_file:
                line_count += 1
        file_list = []
        if line_count/4 > 2000000:
            try:
                os.mkdir(os.path.join(input_folder, genome_name))
            except:
                pass
			#this kind of call forces subprocess to be sequential
			# the real way to do parallel subprocess is subprocess.Popen(<program call>)
            subprocess.call([Python, File_splitter, "2000000", os.path.join(input_folder, genome), os.path.join(input_folder, genome_name)])
            subprocess.call([Python, File_splitter, "2000000", os.path.join(input_folder, genome_name + "2.fastq"), os.path.join(input_folder, genome_name)])
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
        else:
            file_list.append(os.path.join(input_folder, genome_name))
        for item in file_list:
            original = item
            base_name = os.path.splitext(os.path.basename(item))[0]
            
            Input_File = os.path.join(output_folder, base_name, base_name + ".fastq")
            Input_Filepath = os.path.splitext(Input_File)[0]
            Input_File1 = Input_Filepath + "1"
            Input_File2 = Input_Filepath + "2"
            Input_Path = os.path.dirname(Input_File)
            Input_FName = os.path.basename(Input_File)

            try:
                os.chdir(Input_Path)
            except:
                os.mkdir(Input_Path)
                os.chdir(Input_Path)
            
            try:
                os.symlink(original + "1.fastq", os.path.join(output_folder, base_name, base_name + "1.fastq"))
                os.symlink(original + "2.fastq", os.path.join(output_folder, base_name, base_name + "2.fastq"))
            except:
                pass
            
            try:
                File_stats = subprocess.check_output([os.path.join(BBMap_Dir, "testformat.sh"), os.path.join(input_folder, genome_name + "1.fastq")])
            except:
                sys.exit("Use module load java before running script")
            if File_stats.startswith("sanger"):
                Qual = "33"
            else:
                Qual = "64"
            Length = File_stats.split("\t")[4].split("bp")[0]
            
            Host_Contaminants = Input_Filepath + "_host_contaminants_seq.fasta"
            Vector_Contaminants = Input_Filepath + "_vector_contaminants_seq.fasta"
            # Preprocessing
            COMMANDS_Pre = [
            #Python + " " + Sort_Reads + " " + original + "1.fastq" + " " + os.path.join(output_folder, base_name + "1.fastq"),
            #Python + " " + Sort_Reads + " " + original + "2.fastq" + " " + os.path.join(output_folder, base_name + "2.fastq"),
            AdapterRemoval + " --file1 " + Input_File1 + ".fastq" + " --file2 " + Input_File2 + ".fastq" + " --qualitybase " + Qual + " --threads " + Threads + " --minlength " + "30" + " --basename " + os.path.splitext(Input_FName)[0] + "_AdapterRemoval" + " --trimqualities " + " --output1 " + Input_File1 + "_trimmed.fastq" + " --output2 " + Input_File2 + "_trimmed.fastq" + " --singleton " + Input_Filepath + "_singletons_trimmed.fastq",
            vsearch + " --fastq_mergepairs " + Input_File1 + "_trimmed.fastq" + " --reverse " + Input_File2 + "_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastqout " + Input_Filepath + "_overlap_trimmed.fastq" + " --fastqout_notmerged_fwd " + Input_File1 + "_paired_trimmed.fastq" + " --fastqout_notmerged_rev " + Input_File2 + "_paired_trimmed.fastq",
            "cat " + Input_Filepath + "_overlap_trimmed.fastq" + " " + Input_Filepath + "_singletons_trimmed.fastq" + " > " + Input_Filepath + "_unpaired_trimmed.fastq",
            vsearch + " --fastq_filter " + Input_Filepath + "_unpaired_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastq_maxee " + "2.0" + " --fastqout " + Input_Filepath + "_unpaired_quality.fastq",
            vsearch + " --fastq_filter " + Input_File1 + "_paired_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastq_maxee " + "2.0" + " --fastqout " + Input_File1 + "_quality.fastq",
            vsearch + " --fastq_filter " + Input_File2 + "_paired_trimmed.fastq" + " --fastq_ascii " + Qual + " --fastq_maxee " + "2.0" + " --fastqout " + Input_File2 + "_quality.fastq",
            Python + " " + Paired_Reads_Filter + " " + Input_File1 + "_quality.fastq" + " " + Input_File1 + "_paired_quality.fastq" + " " + Input_File2 + "_quality.fastq" + " " + Input_File2 + "_paired_quality.fastq" + " " + Input_Filepath + "_unpaired_quality.fastq",
            cdhit_dup + " -i " + Input_Filepath + "_unpaired_quality.fastq" + " -o " + Input_Filepath + "_unpaired_unique.fastq",
            "mv " + Input_Filepath + "_unpaired_unique.fastq.clstr" + " " + Input_Filepath + "_unpaired.clstr",
            cdhit_dup + " -i " + Input_File1 + "_paired_quality.fastq" + " -i2 " + Input_File2 + "_paired_quality.fastq" + " -o " + Input_File1 + "_paired_unique.fastq" + " -o2 " + Input_File2 + "_paired_unique.fastq",
            "mv " + Input_File1 + "_paired_unique.fastq.clstr" + " " + Input_Filepath + "_paired.clstr",
            "cp " + Host + " " + Host_Contaminants,
            BWA + " index -a bwtsw " + Host_Contaminants,
            SAMTOOLS + " faidx " + Host_Contaminants,
            BWA + " mem -t " + Threads + " " + Host_Contaminants + " " + Input_Filepath + "_unpaired_unique.fastq" + " > " + Input_Filepath + "_unpaired_host_contaminants.sam",
            SAMTOOLS + " view -bS " + Input_Filepath + "_unpaired_host_contaminants.sam" + " > " + Input_Filepath + "_unpaired_host_contaminants.bam",
            SAMTOOLS + " fastq -n -f 4" + " -0 " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_host_contaminants.bam",
            SAMTOOLS + " fastq -n -F 4" + " -0 " + Input_Filepath + "_unpaired_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_host_contaminants.bam",
            BWA + " mem -t " + Threads + " " + Host_Contaminants + " " + Input_File1 + "_paired_unique.fastq" + " " + Input_File2 + "_paired_unique.fastq" + " > " + Input_Filepath + "_paired_host_contaminants.sam",
            SAMTOOLS + " view -bS " + Input_Filepath + "_paired_host_contaminants.sam" + " > " + Input_Filepath + "_paired_host_contaminants.bam",
            SAMTOOLS + " fastq -n -f 13" + " -1 " + Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " -2 " + Input_File2 + "_paired_n_BWA_host_contaminants.fastq" + " " + Input_Filepath + "_paired_host_contaminants.bam",
            SAMTOOLS + " fastq -n -F 4" + " -1 " + Input_File1 + "_paired_host_contaminants.fastq" + " -2 " + Input_File2 + "_paired_host_contaminants.fastq" + " " + Input_Filepath + "_paired_host_contaminants.bam",
            Makeblastdb + " -in " + Host_Contaminants + " -dbtype nucl",
            vsearch + " --fastq_filter " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fasta",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + Host_Contaminants + " " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired_n_host_contaminants.blatout",
            vsearch + " --fastq_filter " + Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_File1 + "_paired_n_BWA_host_contaminants.fasta",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + Host_Contaminants + " " + Input_File1 + "_paired_n_BWA_host_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired_n_host_contaminants.blatout",
            Python + " " + BLAT_Contaminant_Filter + " " + Input_Filepath + "_unpaired_n_BWA_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_n_host_contaminants.blatout" + " " + Input_Filepath + "_unpaired_n_host_contaminants.fastq" + " " + Input_Filepath + "_unpaired_BLAT_host_contaminants.fastq",
            Python + " " + BLAT_Contaminant_Filter + " " + Input_File1 + "_paired_n_BWA_host_contaminants.fastq" + " " + Input_File1 + "_paired_n_host_contaminants.blatout" + " " + Input_File1 + "_paired_n_host_contaminants.fastq" + " " + Input_File1 + "_paired_BLAT_host_contaminants.fastq",
            Python + " " + BLAT_Contaminant_Filter + " " + Input_File2 + "_paired_n_BWA_host_contaminants.fastq" + " " + Input_File1 + "_paired_n_host_contaminants.blatout" + " " + Input_File2 + "_paired_n_host_contaminants.fastq" + " " + Input_File2 + "_paired_BLAT_host_contaminants.fastq",
            "cp " + UniVec_Core + " " + Vector_Contaminants,
            BWA + " index -a bwtsw " + Vector_Contaminants,
            SAMTOOLS + " faidx " + Vector_Contaminants,
            BWA + " mem -t " + Threads + " " + Vector_Contaminants + " " + Input_Filepath + "_unpaired_n_host_contaminants.fastq" + " > " + Input_Filepath + "_unpaired_vector_contaminants.sam",
            SAMTOOLS + " view -bS " + Input_Filepath + "_unpaired_vector_contaminants.sam" + " > " + Input_Filepath + "_unpaired_vector_contaminants.bam",
            SAMTOOLS + " fastq -n -f 4" + " -0 " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " " + Input_Filepath + "_unpaired_vector_contaminants.bam",
            SAMTOOLS + " fastq -n -F 4" + " -0 " + Input_Filepath + "_unpaired_vector_contaminants.fastq" + " " + Input_Filepath + "_unpaired_vector_contaminants.bam",
            BWA + " mem -t " + Threads + " " + Vector_Contaminants + " " + Input_File1 + "_paired_n_host_contaminants.fastq" + " " + Input_File2 + "_paired_n_host_contaminants.fastq" + " > " + Input_Filepath + "_paired_vector_contaminants.sam",
            SAMTOOLS + " view -bS " + Input_Filepath + "_paired_vector_contaminants.sam" + " > " + Input_Filepath + "_paired_vector_contaminants.bam",
            SAMTOOLS + " fastq -n -f 13" + " -1 " + Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " -2 " + Input_File2 + "_paired_n_BWA_vector_contaminants.fastq" + " " + Input_Filepath + "_paired_vector_contaminants.bam",
            SAMTOOLS + " fastq -n -F 4" + " -1 " + Input_File1 + "_paired_vector_contaminants.fastq" + " -2 " + Input_File2 + "_paired_vector_contaminants.fastq" + " " + Input_Filepath + "_paired_vector_contaminants.bam",
            Makeblastdb + " -in " + Vector_Contaminants + " -dbtype nucl",
            vsearch + " --fastq_filter " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fasta",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + Vector_Contaminants + " " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired_n_vector_contaminants.blatout",
            vsearch + " --fastq_filter " + Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " --fastq_ascii " + Qual + " --fastaout " + Input_File1 + "_paired_n_BWA_vector_contaminants.fasta",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + Vector_Contaminants + " " + Input_File1 + "_paired_n_BWA_vector_contaminants.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired_n_vector_contaminants.blatout",
            Python + " " + BLAT_Contaminant_Filter + " " + Input_Filepath + "_unpaired_n_BWA_vector_contaminants.fastq" + " " + Input_Filepath + "_unpaired_n_vector_contaminants.blatout" + " " + Input_Filepath + "_unpaired_n_contaminants.fastq" + " " + Input_Filepath + "_unpaired_BLAT_vector_contaminants.fastq",
            Python + " " + BLAT_Contaminant_Filter + " " + Input_File1 + "_paired_n_BWA_vector_contaminants.fastq" + " " + Input_File1 + "_paired_n_vector_contaminants.blatout" + " " + Input_File1 + "_paired_n_contaminants.fastq" + " " + Input_File1 + "_paired_BLAT_vector_contaminants.fastq",
            Python + " " + BLAT_Contaminant_Filter + " " + Input_File2 + "_paired_n_BWA_vector_contaminants.fastq" + " " + Input_File1 + "_paired_n_vector_contaminants.blatout" + " " + Input_File2 + "_paired_n_contaminants.fastq" + " " + Input_File2 + "_paired_BLAT_vector_contaminants.fastq",
            "mkdir -p " + os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants",
            Python + " " + File_splitter + " " + "10000" + " " + Input_Filepath + "_unpaired_n_contaminants.fastq" + " " + os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants",
            "mkdir -p " + os.path.splitext(Input_FName)[0] + "_paired_n_contaminants",
            Python + " " + File_splitter + " " + "10000" + " " + Input_File1 + "_paired_n_contaminants.fastq" + " " + os.path.splitext(Input_FName)[0] + "_paired_n_contaminants",
            Python + " " + File_splitter + " " + "10000" + " " + Input_File2 + "_paired_n_contaminants.fastq" + " " + os.path.splitext(Input_FName)[0] + "_paired_n_contaminants"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Preprocess.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Preprocess")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Pre))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Pre = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Preprocess.pbs"])

            COMMANDS_rRNA = [
            "JOBS=$(" + Python + " " + rRNA_Split_Jobs + " " + Input_File + " " + JobID_Pre.strip("\n") + ");" + "qalter -W depend=afterok:$JOBS $JOB2"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_rRNA_Submit.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_rRNA_Submit")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_rRNA))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_rRNA = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_rRNA_Submit.pbs", "-W", "depend=afterok:" + JobID_Pre.strip("\n")])

            COMMANDS_Combine = [
            "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.basename(Input_Filepath) + "_unpaired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + Input_Filepath + "_mRNA_unpaired.fastq",
            "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_unpaired_n_contaminants", os.path.basename(Input_Filepath) + "_unpaired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + Input_Filepath + "_rRNA_unpaired.fastq",
            "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File1) + "_paired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + Input_File1 + "_mRNA.fastq",
            "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File1) + "_paired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + Input_File1 + "_rRNA.fastq",
            "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File2) + "_paired_n_contaminants" + "_split_*" + "_mRNA.fastq") + " > " + Input_File2 + "_mRNA.fastq",
            "cat " + os.path.join(os.path.splitext(Input_FName)[0] + "_paired_n_contaminants", os.path.basename(Input_File2) + "_paired_n_contaminants" + "_split_*" + "_rRNA.fastq") + " > " + Input_File2 + "_rRNA.fastq",
            Python + " " + Reduplicate + " " + Input_Filepath + "_unpaired_quality.fastq" + " " + Input_Filepath + "_mRNA_unpaired.fastq" + " " + Input_Filepath + "_unpaired.clstr" + " " + Input_Filepath + "_all_mRNA_unpaired.fastq",
            Python + " " + Reduplicate + " " + Input_File1 + "_paired_quality.fastq" + " " + Input_File1 + "_mRNA.fastq" + " " + Input_Filepath + "_paired.clstr" + " " + Input_File1 + "_all_mRNA.fastq",
            Python + " " + Reduplicate + " " + Input_File2 + "_paired_quality.fastq" + " " + Input_File2 + "_mRNA.fastq" + " " + Input_Filepath + "_paired.clstr" + " " + Input_File2 + "_all_mRNA.fastq"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Combine.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Combine")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Combine))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Combine = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Combine.pbs", "-W", "depend=afterok:" + JobID_rRNA.strip("\n")])

            subprocess.call(["qalter", "-v", "JOB2=" + JobID_Combine.strip("\n").split(".")[0], JobID_rRNA.strip("\n")])

            #BIGG Database, AGORA Nature paper, Additional functionality
            # Transcript Assembly
            Contigs = os.path.join(Input_Path, os.path.splitext(Input_FName)[0] + "_SpadesOut", "contigs.fasta")
            COMMANDS_Assemble = [
            "mkdir -p " + os.path.join(Input_Path, os.path.splitext(Input_FName)[0]) + "_SpadesOut",
            #Python + " " + Spades + " --rna" + " -1 " + Input_File1 + "_all_mRNA.fastq" + " -2 " + Input_File2 + "_all_mRNA.fastq" + " -s " + Input_Filepath + "_all_mRNA_unpaired.fastq" + " -o " + os.path.splitext(Input_FName)[0] + "_SpadesOut",
            Python + " " + Spades + " -k 21,33,55,77 --meta" + " -1 " + Input_File1 + "_all_mRNA.fastq" + " -2 " + Input_File2 + "_all_mRNA.fastq" + " -o " + os.path.splitext(Input_FName)[0] + "_SpadesOut",
            BWA + " index -a bwtsw " + Contigs,
            BWA + " mem -t " + Threads + " -B 40 -O 60 -E 10 -L 50 " + Contigs + " " + Input_File1 + "_all_mRNA.fastq " + Input_File2 + "_all_mRNA.fastq | " + SAMTOOLS + " view > " + Input_Filepath + "_contig_paired.sam",
            BWA + " mem -t " + Threads + " -B 40 -O 60 -E 10 -L 50 " + Contigs + " " + Input_Filepath + "_all_mRNA_unpaired.fastq | " + SAMTOOLS + " view > " + Input_Filepath + "_contig_unpaired.sam",
            Python + " " + Map_reads_contigs + " " + Input_File1 + "_all_mRNA.fastq" + " " + Input_File2 + "_all_mRNA.fastq" + " " + Input_Filepath + "_all_mRNA_unpaired.fastq" + " " + Input_Filepath + "_contig_paired.sam" + " " + Input_Filepath + "_contig_unpaired.sam" + " " + Input_Filepath + "_contig_map.tsv"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Assemble.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Assemble")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Assemble))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Assemble = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Assemble.pbs", "-W", "depend=afterok:" + JobID_Combine.strip("\n")])

            # Protein Annotation BWA
            COMMANDS_Annotate_BWA = [
            BWA + " mem -t " + Threads + " " + DNA_DB + " " + Contigs + " | " + SAMTOOLS + " view > " + Input_Filepath + "_contigs_BWA.sam",
            BWA + " mem -t " + Threads + " " + DNA_DB + " " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " | " + SAMTOOLS + " view > " + Input_Filepath + "_unpaired_unmapped_BWA.sam",
            BWA + " mem -t " + Threads + " " + DNA_DB + " " + Input_File1 + "_all_mRNA_unmapped.fastq" + " " + Input_File2 + "_all_mRNA_unmapped.fastq" + " | " + SAMTOOLS + " view > " + Input_Filepath + "_paired_unmapped_BWA.sam",
            Python + " " + Map_reads_gene_BWA + " " + DNA_DB + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Contigs + " " + Input_Filepath + "_contigs_BWA.sam" + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " " + Input_Filepath + "_unpaired_unmapped_BWA.sam" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " " + Input_File1 + "_all_mRNA_unmapped.fastq" + " " + Input_Filepath + "_paired_unmapped_BWA.sam" + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " " + Input_File2 + "_all_mRNA_unmapped.fastq" + " " + Input_Filepath + "_paired_unmapped_BWA.sam" + " " + Input_File2 + "_unmapped_n_BWA.fasta",
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_BWA.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_HighMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_BWA")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_BWA))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_BWA = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BWA.pbs", "-W", "depend=afterok:" + JobID_Assemble.strip("\n")])

            # Protein Annotation BLAT 1
            COMMANDS_Annotate_BLAT1 = [
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_1" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_1" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_2" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_2" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_3" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_3" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_4" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_4" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_5" + DNA_DB_Extension + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_contigs" + "_5" + ".blatout",
            "cat " + Input_Filepath + "_contigs" + "_[1-5]" + ".blatout" + " > " + Input_Filepath + "_contigs" + ".blatout"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_BLAT1.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_HighMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT1")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_BLAT1))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_BLAT1 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT1.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT 2
            COMMANDS_Annotate_BLAT2 = [
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_1" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_1" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_2" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_2" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_3" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_3" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_4" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_4" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_5" + DNA_DB_Extension + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_Filepath + "_unpaired" + "_5" + ".blatout",
            "cat " + Input_Filepath + "_unpaired" + "_[1-5]" + ".blatout" + " > " + Input_Filepath + "_unpaired" + ".blatout"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_BLAT2.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_HighMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT2")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_BLAT2))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_BLAT2 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT2.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT 3
            COMMANDS_Annotate_BLAT3 = [
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_1" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_1" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_2" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_2" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_3" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_3" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_4" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_4" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_5" + DNA_DB_Extension + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File1 + "_paired" + "_5" + ".blatout",
            "cat " + Input_File1 + "_paired" + "_[1-5]" + ".blatout" + " > " + Input_File1 + "_paired" + ".blatout"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_BLAT3.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_HighMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT3")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_BLAT3))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_BLAT3 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT3.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT 4
            COMMANDS_Annotate_BLAT4 = [
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_1" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_1" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_2" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_2" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_3" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_3" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_4" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_4" + ".blatout",
            BLAT + " -noHead -minIdentity=90 -minScore=65 " + DNA_DB_Prefix + "_5" + DNA_DB_Extension + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " -fine -q=rna -t=dna -out=blast8 -threads=" + Threads + " " + Input_File2 + "_paired" + "_5" + ".blatout",
            "cat " + Input_File2 + "_paired" + "_[1-5]" + ".blatout" + " > " + Input_File2 + "_paired" + ".blatout"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_BLAT4.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_HighMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT4")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_BLAT4))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_BLAT4 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT4.pbs", "-W", "depend=afterok:" + JobID_Annotate_BWA.strip("\n")])

            # Protein Annotation BLAT Postprocessing
            COMMANDS_Annotate_BLAT_Post = [
            Python + " " + Map_reads_gene_BLAT + " " + DNA_DB + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_genes.fna" + " " + Input_Filepath + "_contigs_n_BWA.fasta" + " " + Input_Filepath + "_contigs" + ".blatout" + " " + Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA.fasta" + " " + Input_Filepath + "_unpaired" + ".blatout" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " " + Input_File1 + "_unmapped_n_BWA.fasta" + " " + Input_File1 + "_paired" + ".blatout" + " " + Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " " + Input_File2 + "_unmapped_n_BWA.fasta" + " " + Input_File2 + "_paired" + ".blatout" + " " + Input_File2 + "_unmapped_n_BWA_BLAT.fasta"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_BLAT_Postprocessing.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT_Postprocessing")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_BLAT_Post))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_BLAT_Post = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_BLAT_Postprocessing.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT1.strip("\n") + ":" + JobID_Annotate_BLAT2.strip("\n") + ":" + JobID_Annotate_BLAT3.strip("\n") + ":" + JobID_Annotate_BLAT4.strip("\n")])


            # Protein Annotation Diamond 1
            COMMANDS_Annotate_Diamond1 = [
            "mkdir -p " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp1",
            DIAMOND + " blastx -p " + Threads + " -d " + Prot_DB + " -q " + Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " -o " + Input_Filepath + "_contigs.dmdout" + " -f 6 -t " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp1 -k 10 --id 85 --query-cover 65 --min-score 60",
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_DMD1.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_DMD1")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_Diamond1))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_Diamond1 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD1.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond 2
            COMMANDS_Annotate_Diamond2 = [
            "mkdir -p " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp2",
            DIAMOND + " blastx -p " + Threads + " -d " + Prot_DB + " -q " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " -o " + Input_Filepath + "_unpaired.dmdout" + " -f 6 -t " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp2 -k 10 --id 85 --query-cover 65 --min-score 60",
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_DMD2.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_DMD2")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_Diamond2))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_Diamond2 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD2.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond 3
            COMMANDS_Annotate_Diamond3 = [
            "mkdir -p " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp3",
            DIAMOND + " blastx -p " + Threads + " -d " + Prot_DB + " -q " + Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " -o " + Input_File1 + "_paired.dmdout" + " -f 6 -t " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp3 -k 10 --id 85 --query-cover 65 --min-score 60",
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_DMD3.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_DMD3")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_Diamond3))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_Diamond3 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD3.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond 4
            COMMANDS_Annotate_Diamond4 = [
            "mkdir -p " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp4",
            DIAMOND + " blastx -p " + Threads + " -d " + Prot_DB + " -q " + Input_File2 + "_unmapped_n_BWA_BLAT.fasta" + " -o " + Input_File2 + "_paired.dmdout" + " -f 6 -t " + os.path.splitext(Input_FName)[0] + "_dmnd_tmp4 -k 10 --id 85 --query-cover 65 --min-score 60"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_DMD4.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_DMD4")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_Diamond4))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_Diamond4 = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD4.pbs", "-W", "depend=afterok:" + JobID_Annotate_BLAT_Post.strip("\n")])

            # Protein Annotation Diamond Postprocess
            COMMANDS_Annotate_Diamond_Post = [
            Python + " " + Map_reads_prot_DMND + " " + Prot_DB + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_genes.fna" + " " + Input_Filepath + "_proteins.faa" + " " + Input_Filepath + "_contigs_n_BWA_BLAT.fasta" + " " + Input_Filepath + "_contigs.dmdout" + " " + Input_Filepath + "_contigs_n_BWA_BLAT_DMD.fasta" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT.fasta" + " " + Input_Filepath + "_unpaired.dmdout" + " " + Input_Filepath + "_unpaired_unmapped_n_BWA_BLAT_DMD.fasta" + " " + Input_File1 + "_unmapped_n_BWA_BLAT.fasta" + " " + Input_File1 + "_paired.dmdout" + " " + Input_File1 + "_unmapped_n_BWA_BLAT_DMD.fasta" + " " + Input_File2 + "_unmapped_n_BWA_BLAT.fasta" + " " + Input_File2 + "_paired.dmdout" + " " + Input_File2 + "_unmapped_n_BWA_BLAT_DMD.fasta"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Annotate_DMD_Postprocess.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Annotate_DMD_Postprocess")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Annotate_Diamond_Post))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Annotate_Diamond_Post = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Annotate_DMD_Postprocess.pbs", "-W", "depend=afterany:" + JobID_Annotate_Diamond1.strip("\n") + ":" + JobID_Annotate_Diamond2.strip("\n") + ":" + JobID_Annotate_Diamond3.strip("\n") + ":" + JobID_Annotate_Diamond4.strip("\n")])
            
            # Classify Reads
            COMMANDS_Classify = [
            Python + " " + Annotated_taxid + " " + Input_Filepath + "_gene_map.tsv" + " " + accession2taxid + " " + Input_Filepath + "_TaxIDOut.tsv",
            Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + Contigs + " -z " + Threads + " -o " + Input_Filepath + "_contigs_KaijuOut.tsv",
            Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " -z " + Threads + " -o " + Input_Filepath + "_unpaired_KaijuOut.tsv",
            Kaiju + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -f " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/kaiju_db_nr.fmi" + " -i " + Input_File1 + "_all_mRNA_unmapped.fastq" + " -j " + Input_File2 + "_all_mRNA_unmapped.fastq" + " -z " + Threads + " -o " + Input_Filepath + "_paired_KaijuOut.tsv",
            "cat " + Input_Filepath + "_contigs_KaijuOut.tsv" + " " + Input_Filepath + "_unpaired_KaijuOut.tsv" + " " + Input_Filepath + "_paired_KaijuOut.tsv" + " > " + Input_Filepath + "_KaijuOut.tsv",
            Centrifuge + " -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -1 " + Input_File1 + "_all_mRNA_unmapped.fastq" + " -2 " + Input_File2 + "_all_mRNA_unmapped.fastq" + " -U " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq" + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + Qual + " -p " + Threads + " -S " + Input_Filepath + "_unmapped_CentrifugeOut.tsv" + " --report-file " + Input_Filepath + "_unmapped_CentrifugeReport.txt",
            Centrifuge + " -f -x " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nt" + " -U " + Contigs + " --exclude-taxids 2759 --tab-fmt-cols " + "score,readID,taxID" + " --phred" + Qual + " -p " + Threads + " -S " + Input_Filepath + "_contigs_CentrifugeOut.tsv" + " --report-file " + Input_Filepath + "_contigs_CentrifugeReport.txt",
            "cat " + Input_Filepath + "_unmapped_CentrifugeOut.tsv" + " " + Input_Filepath + "_contigs_CentrifugeOut.tsv" + " > " + Input_Filepath + "_CentrifugeOut.tsv",
            #kSLAM + " " + "--db=/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/" + " " + "--output-file=" + Input_Filepath + "_unpaired_kSLAMOut" + " " + Input_Filepath + "_all_mRNA_unpaired_unmapped.fastq",
            #kSLAM + " " + "--db=/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/" + " " + "--output-file=" + Input_Filepath + "_paired_kSLAMOut" + " " + Input_File1 + "_all_mRNA_unmapped.fastq" + " " + Input_File2 + "_all_mRNA_unmapped.fastq",
            #"sed \'s/^.*/C\\t&/\' " + Input_Filepath + "_unpaired_kSLAMOut_PerRead" + " " + Input_Filepath + "_paired_kSLAMOut_PerRead" + " > " + Input_Filepath + "_kSLAMOut.tsv",
            Python + " " + Classification_combine + " " + Input_Filepath + "_contig_map.tsv" + " " + Input_Filepath + "_WEVOTEOut_ensemble.csv" + " " + Input_Filepath + "_TaxIDOut.tsv" + " " + Input_Filepath + "_KaijuOut.tsv" + " " + Input_Filepath + "_CentrifugeOut.tsv",# + " " + Input_Filepath + "_kSLAMOut.tsv",
            "mkdir -p " + Input_Filepath + "_WEVOTEOut",
            "cp "  + Input_Filepath + "_WEVOTEOut_ensemble.csv" + " " + Input_Filepath + "_WEVOTEOut",
            "cd " + os.path.dirname(WEVOTE),
            WEVOTE + " -o " + Input_Filepath + "_WEVOTEOut" + " --db " + WEVOTEDB + " -c",
            "cd $PBS_O_WORKDIR",
            "awk -F \'\\t\' \'{print \"C\\t\"$1\"\\t\"$9}\' " + os.path.join(Input_Filepath + "_WEVOTEOut", os.path.splitext(Input_FName)[0] + "_WEVOTEOut_WEVOTE_Details.txt") + " > " + Input_Filepath + "_WEVOTEOut.tsv",
            Python + " " + Contrain_classification + " " + "family" + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Nodes + " " + Names + " " + Input_Filepath + "_WEVOTEOut_family.tsv",
            Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + Input_Filepath + "_WEVOTEOut_family.tsv" + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.txt",
            "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
            ktImportText + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Classify.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_vHighMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Classify")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Classify))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Classify = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Classify.pbs", "-W", "depend=afterok:" + JobID_Annotate_Diamond_Post.strip("\n")])
            
            # Prepare EC annotation files
            EC_Split = os.path.join(Input_Filepath + "_EC_Annotation", "Split")
            EC_Output = os.path.join(Input_Filepath + "_EC_Annotation", "Output")
            COMMANDS_EC_Preprocess = [
            "mkdir -p " + EC_Split,
            "mkdir -p " + EC_Output,
            Python + " " + File_splitter + " " + "1000" + " " + Input_Filepath + "_proteins.faa" + " " + EC_Split,
            ]

            with open(os.path.splitext(Input_FName)[0] + "_EC_Preprocess.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_EC_Preprocess")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_EC_Preprocess))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_EC_Preprocess = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_EC_Preprocess.pbs", "-W", "depend=afterok:" + JobID_Annotate_Diamond_Post.strip("\n")])

            COMMANDS_Detect = [
            "JOBS=$(" + Python + " " + Detect_Submit + " " + EC_Split + " " + EC_Output + " " + Threads + " " + JobID_EC_Preprocess.strip("\n") + ");" + "qalter -W depend=afterok:$JOBS $JOB2"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Detect.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Detect")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Detect))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Detect = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Detect.pbs", "-W", "depend=afterok:" + JobID_EC_Preprocess.strip("\n")])

            COMMANDS_Combine_Detect = ["cat " + os.path.join(EC_Output, "Detect", "*.toppred") + " > " + os.path.join(EC_Output, "Detect", os.path.splitext(Input_FName)[0] + "_proteins.toppred")]

            with open(os.path.splitext(Input_FName)[0] + "_Combine_Detect.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Combine_Detect")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Combine_Detect))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Combine_Detect = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Combine_Detect.pbs", "-W", "depend=afterok:" + JobID_Detect.strip("\n")])

            subprocess.call(["qalter", "-v", "JOB2=" + JobID_Combine_Detect.strip("\n").split(".")[0], JobID_Detect.strip("\n")])

            COMMANDS_PRIAM = [
            "mkdir -p " + os.path.join(EC_Output, "PRIAM"),
            "cd " + os.path.join(EC_Output, "PRIAM"),
            "java -jar" + " " + Priam + " -n " + os.path.splitext(Input_FName)[0] + "_PRIAM" + " -i " + Input_Filepath + "_proteins.faa" + " -p " + os.path.join(os.path.dirname(Priam), "PRIAM_MAR15") + " -od " + os.path.join(EC_Output, "PRIAM") +" -e T -pt 0.5 -mo -1 -mp 70 -cc T -cg T -bd " + BLAST_dir,
            ]

            with open(os.path.splitext(Input_FName)[0] + "_PRIAM.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_PRIAM")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_PRIAM))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_PRIAM = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_PRIAM.pbs", "-W", "depend=afterok:" + JobID_EC_Preprocess.strip("\n")])

            COMMANDS_EC_Diamond = [
            "mkdir -p " + os.path.join(EC_Output, "Diamond"),
            "cd " + os.path.join(EC_Output, "Diamond"),
            DIAMOND + " blastp -p " + Threads + " --query "+ Input_Filepath + "_proteins.faa" + " --db "+ SWISS_PROT + " --outfmt "+ "6 qseqid sseqid qstart qend sstart send evalue bitscore qcovhsp slen pident" + " --out " + os.path.join(EC_Output, "Diamond", os.path.splitext(Input_FName)[0] + ".blastout") + " --evalue 0.0000000001 --max-target-seqs 1"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_EC_Diamond.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_EC_Diamond")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_EC_Diamond))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_EC_Diamond = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_EC_Diamond.pbs", "-W", "depend=afterok:" + JobID_EC_Preprocess.strip("\n")])

            # EC Annotation Compile
            COMMANDS_EC_Postprocess = [
            Python + " " + EC_Annotation_Post + " " + Input_Filepath + "_proteins.faa" + " " + EC_Output
            #Produce_Table
            ]

            with open(os.path.splitext(Input_FName)[0] + "_EC_Postprocess.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_EC_Postprocess")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_EC_Postprocess))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_EC_Postprocess = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_EC_Postprocess.pbs", "-W", "depend=afterok:" + JobID_Combine_Detect.strip("\n") + ":" + JobID_PRIAM.strip("\n") + ":" + JobID_EC_Diamond.strip("\n")])

            # Network Generation
            COMMANDS_Network = [
            Python + " " + RPKM + " " + Nodes + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + os.path.join(Input_Filepath + "_EC_Annotation", "Output", "Consolidated", os.path.splitext(Input_FName)[0] + "_proteins.ECs_All") + " " + Input_Filepath + "_RPKM.tsv" + " " + Input_Filepath + "_Cytoscape.tsv"
            ]

            with open(os.path.splitext(Input_FName)[0] + "_Network.pbs", "w") as PBS_script_out:
                for line in PBS_Submit_LowMem.splitlines():
                    if "NAME" in line:
                        line = line.replace("NAME", os.path.splitext(Input_FName)[0] + "_Network")
                    if "COMMANDS" in line:
                        PBS_script_out.write("\n".join(COMMANDS_Network))
                        break
                    PBS_script_out.write(line + "\n")
            JobID_Network = subprocess.check_output(["qsub", os.path.splitext(Input_FName)[0] + "_Network.pbs", "-W", "depend=afterok:" + JobID_EC_Postprocess.strip("\n") + ":" + JobID_Classify.strip("\n")])
            Network_list.append(JobID_Network.strip("\n"))

    if len(file_list) > 1:
        os.chdir(output_folder)
        Input_Filepath = os.path.join(output_folder, os.path.splitext(os.path.basename(file_list[0]))[0])[:-11]
        COMMANDS_Join = [
        "cat " + Input_Filepath + "*/*_WEVOTEOut.tsv" + " > " + Input_Filepath + "_WEVOTEOut.tsv",
        Python + " " + Contrain_classification + " " + "family" + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Nodes + " " + Names + " " + Input_Filepath + "_WEVOTEOut_family.tsv",
        Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + Input_Filepath + "_WEVOTEOut_family.tsv" + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.txt",
        "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
        ktImportText + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
        "cat " + Input_Filepath + "*/*_gene_map.tsv" + " > " + Input_Filepath + "_gene_map.tsv",
        "cat " + os.path.join(Input_Filepath + "*/*_EC_Annotation", "Output", "Consolidated", "*_proteins.ECs_All") + " > " + Input_Filepath + "_proteins.ECs_All",
        Python + " " + RPKM + " " + Nodes + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_proteins.ECs_All" + " " + Input_Filepath + "_RPKM.tsv" + " " + Input_Filepath + "_Cytoscape.tsv"
        ]

        with open(os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join.pbs", "w") as PBS_script_out:
            for line in PBS_Submit_LowMem.splitlines():
                if "NAME" in line:
                    line = line.replace("NAME", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join")
                if "COMMANDS" in line:
                    PBS_script_out.write("\n".join(COMMANDS_Join))
                    break
                PBS_script_out.write(line + "\n")
        JobID_Join = subprocess.check_output(["qsub", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join.pbs", "-W", "depend=afterok:" + ":".join(Network_list)])