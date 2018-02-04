#!/bin/bash
# Script for running mt_joinsplits_RPKM.py

# "output" folder must have subfolders containing splitfiles

module purge
module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras python

OUTPUT_FOLDER=$SCRATCH/datasets/Ilott_2016/LARGEout
JOIN_SCRIPT=$HOME/parkinson_int_meta_pipe/jordans_scripts/python/mt_joinsplits_RPKM.py

python $JOIN_SCRIPT $OUTPUT_FOLDER
