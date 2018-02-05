#!/bin/bash
# Script for running mt_splitfiles.py

# This script will stay active until ALL the .fasta files in the input folder are split up.

module purge
module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras python

INPUT_FOLDER=$SCRATCH/datasets/Ilott_2016/50Kin
SPLIT_SCRIPT=$HOME/parkinson_int_meta_pipe/jordans_scripts/python/mt_splitfiles.py
SPLIT_READS_NUM=1000000

python $SPLIT_SCRIPT $INPUT_FOLDER $SPLIT_READS_NUM
