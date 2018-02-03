#!/bin/bash
# Script for running mt_pipe_isolate.py

# Input folder must have subfolders containing splitfiles

# All jobs get submitted into the qsub queue on the fly, but the script will close before the jobs complete. The jobs have built-in dependency for start times.

# Job order (brackets not start points):
# 'Pre', 'Isolate', 'rRNA', ('Combine'), 'Assemble', 'Annotate_BWA', 'Annotate_BWA_Post', 'Annotate_BLAT', 'Annotate_BLAT_Post', 'Annotate_Diamond', 'Annotate_Diamond_Post', 'Classify', 'EC_Preprocess', 'Detect', ('Combine_Detect'), 'PRIAM', 'EC_Diamond', 'EC_Postprocess', 'Network']

# START_POINT & END_POINT are inclusive

module purge
module load gcc/5.2.0 boost/1.60.0-gcc5.2.0 intel/15.0.2 openmpi java blast extras python

INPUT_FOLDER=$SCRATCH/datasets/Ilott_2016/100KinTEST
OUTPUT_FOLDER=$SCRATCH/datasets/Ilott_2016/100KoutTEST
PIPELINE_SCRIPT=$HOME/scripts/python/mt_pipe_TEST.py

START_POINT=Pre
END_POINT=Annotate_BLAT

mkdir -p $OUTPUT_FOLDER

python $PIPELINE_SCRIPT $INPUT_FOLDER $OUTPUT_FOLDER $START_POINT $END_POINT
