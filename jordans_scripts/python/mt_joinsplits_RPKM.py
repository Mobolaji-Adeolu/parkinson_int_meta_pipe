#!/usr/bin/env python

# Joining from all split folders in output folder to files in output folder.

import sys
import os
import os.path
import subprocess
import multiprocessing
import glob

cdhit_dup = "/home/j/jparkins/mobolaji/Tools/CDHIT/cd-hit-v4.6.6-2016-0711/cd-hit-auxtools/cd-hit-dup"
Timmomatic = "/home/j/jparkins/mobolaji/Tools/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar"
AdapterRemoval = "/home/j/jparkins/mobolaji/Tools/AdapterRemoval/adapterremoval-2.2.0/build/AdapterRemoval"
vsearch = "/home/j/jparkins/mobolaji/Tools/vsearch/vsearch-2.4.2-linux-x86_64/bin/vsearch"
UniVec_Core = "/home/j/jparkins/mobolaji/Databases/UniVec_Core.fasta"
Flash = "/home/j/jparkins/mobolaji/Tools/Flash/FLASH-1.2.11/flash"
Perl = "/home/j/jparkins/mobolaji/perl"
Perl_Script_Dir = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Xuejian"
#Python = "/home/j/jparkins/mobolaji/python"
Python = "/scinet/gpc/tools/Python/Python272-shared/bin/python"
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
#Map_reads_contigs = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Map_read_contigs.py"
Map_reads_contigs = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/map_read_contigs.py"
Paired_Reads_Filter = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Paired_Reads_Filter.py"
#BLAT_Contaminant_Filter = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/BLAT_Contaminant_Filter.py"
BLAT_Contaminant_Filter = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/BLAT_contaminant_filter.py"
#File_splitter = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/File_splitter.py"
File_splitter = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/file_splitter.py"
Sort_Reads = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/mt_sortreads.py"
# Sort_Reads = "/home/j/jparkins/mobolaji/Read_Classification/Sort_Reads.py"
#rRNA_Split_Jobs = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/rRNA_Split_Jobs.py"
rRNA_Split_Jobs = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/rRNA_split_jobs.py"
#Map_reads_gene_BWA = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Map_read_gene_BWA.py"
Map_reads_gene_BWA = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/map_read_gene_BWA.py"
#Map_reads_gene_BLAT = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Map_read_gene_BLAT.py"
Map_reads_gene_BLAT = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/map_read_gene_BLAT.py"
#Map_reads_prot_DMND = "/home/j/jparkins/mobolaji/Metatranscriptome_Scripts/Mobolaji/Map_read_prot_DMND.py"
Map_reads_prot_DMND = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/map_read_prot_DMND.py"
Spades = "/home/j/jparkins/mobolaji/Tools/SPAdes/SPAdes-3.9.1-Linux/bin/spades.py"
EC_Annotation_Prep = "/home/j/jparkins/mobolaji/EC_Prediction_Scripts/0_Preprocess_Input.py"
#Detect_Submit = "/home/j/jparkins/mobolaji/EC_Prediction_Scripts/1-1a_Detect_Submission.py"
Detect_Submit = "/home/j/jparkins/ang/parkinson_int_meta_pipe/jordans_scripts/python/detect_submission.py"
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

PBS_Submit = """#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l WALLTIME
#PBS -N NAME

module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras python
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8
OLDPATH=$PATH:/home/j/jparkins/ctorma/emboss/bin/:/home/j/jparkins/mobolaji/Tools/Barrnap/bin/:/home/j/jparkins/mobolaji/Tools/HMMer/hmmer-3.1b2-linux-intel-x86_64/binaries/:/home/j/jparkins/mobolaji/Tools/Python27/Python-2.7.12/:/home/j/jparkins/mobolaji/Tools/Bowtie2/bowtie2-2.3.0/:/home/j/jparkins/mobolaji/Tools/SAMTOOLS/samtools-1.3.1/
NEWPATH=/home/j/jparkins/mobolaji/:$OLDPATH
export PATH=$NEWPATH

COMMANDS"""


output_folder = sys.argv[1]

os.chdir(output_folder)

# find all split folders:
split_folders = glob.glob('*_split_*_')

# extract list of genomes with split folders:
genome_join = []
for name in split_folders:
    genome_name = name.split("_split_")[0]
    if genome_name not in genome_join:
        genome_join.append(genome_name)

# do Join.pbs scripts already exist for these genomes?:
pbs_exists = []
for name in genome_join:
    filepath = os.path.join(output_folder, name + "_Join.pbs")
    if os.path.exists(filepath):
        pbs_exists.append(name)

# if scripts exist, list and quit:
if len(pbs_exists) > 0:
    print "Join.pbs scrips already exist for the following genomes:"
    for name in pbs_exists:
        print name[:-1]
    sys.exit()

print 'Joining splits for:'

# create Join.pbs scripts:
for input_file in genome_join:
    Input_Filepath = os.path.join(output_folder, input_file)

    print str(input_file) + ' ...started'

    # join commands:
    COMMANDS_Join = [
    # JOIN WEVOTE TAXA OUTPUTS:
    "cat " + Input_Filepath + "*/*_WEVOTEOut.tsv" + " > " + Input_Filepath + "_WEVOTEOut.tsv",
#    # RECLASSIFY TO FAMILY LEVEL DEPTH:
#    Python + " " + Contrain_classification + " " + "family" + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Nodes + " " + Names + " " + Input_Filepath + "_WEVOTEOut_family.tsv",
#    # GENERATE HIERARCHICAL MULTI-LAYER PIE CHART OF FAMILY COMPOSITION:
#    Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + Input_Filepath + "_WEVOTEOut_family.tsv" + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.txt",
#    "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + Input_Filepath + "_WEVOTEOut_family_Krona.txt" + " > " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
#    ktImportText + " -o " + Input_Filepath + "_WEVOTEOut_family_Krona.html" + " " + Input_Filepath + "_WEVOTEOut_family_Krona.tsv",
    # RECLASSIFY TO SPECIES LEVEL DEPTH:
    Python + " " + Contrain_classification + " " + "species" + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Nodes + " " + Names + " " + Input_Filepath + "_WEVOTEOut_species.tsv",
    # GENERATE HIERARCHICAL MULTI-LAYER PIE CHART OF SPECIES COMPOSITION:
    Kaiju2krona + " -t " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/nodes_nr.dmp" + " -n " + "/scratch/j/jparkins/mobolaji/NCBI_nr_db/Index/names_nr.dmp" + " -i " + Input_Filepath + "_WEVOTEOut_species.tsv" + " -o " + Input_Filepath + "_WEVOTEOut_species_Krona.txt",
    "awk -F \'\\t\' \'{OFS=\"\\t\";$2=\"\";$3=\"\";print}\' " + Input_Filepath + "_WEVOTEOut_species_Krona.txt" + " > " + Input_Filepath + "_WEVOTEOut_species_Krona.tsv",
    ktImportText + " -o " + Input_Filepath + "_WEVOTEOut_species_Krona.html" + " " + Input_Filepath + "_WEVOTEOut_species_Krona.tsv",
    # JOIN ANNOTATED GENES/PROTEINS:
    "cat " + Input_Filepath + "*/*_gene_map.tsv" + " > " + Input_Filepath + "_gene_map.tsv",
    # JOIN ANNOTATED ECs:
    "cat " + os.path.join(Input_Filepath + "*/*_EC_Annotation", "Output", "Consolidated", "*_proteins.ECs_All") + " > " + Input_Filepath + "_proteins.ECs_All",
    # CALCULATE RPKM & COMPILE INFO '_RPKM.tsv'= [geneID/proteinID, length, #reads, EC#, Total RPKM, RPKM per phylum]:
    # COMPILE EC NODE ATTRIBUTES FOR CYTOSCAPE '_Cytoscape.tsv'= [EC#, RPKM, RPKM per phylum]:
    Python + " " + RPKM + " " + Nodes + " " + Input_Filepath + "_WEVOTEOut.tsv" + " " + Input_Filepath + "_gene_map.tsv" + " " + Input_Filepath + "_proteins.ECs_All" + " " + Input_Filepath + "_RPKM.tsv" + " " + Input_Filepath + "_Cytoscape.tsv"
    ]

    # write script:
    with open(os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join.pbs", "w") as PBS_script_out:
        for line in PBS_Submit.splitlines():
            if "WALLTIME" in line:
                        line = line.replace("WALLTIME", "walltime=12:00:00")
            if "NAME" in line:
                line = line.replace("NAME", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join")
            if "COMMANDS" in line:
                PBS_script_out.write("\n".join(COMMANDS_Join))
                break
            PBS_script_out.write(line + "\n")

    # run script:
    JobID_Join = subprocess.check_output(["qsub", os.path.splitext(os.path.basename(Input_Filepath))[0] + "_Join.pbs"])
