#!/bin/bash
# Script for running mt_splitfiles.py

module purge
module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras python

SORT_SCRIPT=$HOME/parkinson_int_meta_pipe/jordans_scripts/python/mt_sortreads.py
UNSORTED=$SCRATCH/datasets/NOD/unsorted/NOD504CecQY_2.fastq
SORTED_FOLDER=$SCRATCH/datasets/NOD/sorted
SORTED=$SORTED_FOLDER/NOD504CecQY_2.fastq

mkdir -p $SORTED_FOLDER

python $SORT_SCRIPT $UNSORTED $SORTED
