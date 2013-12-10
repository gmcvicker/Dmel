#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 1-4
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 2GB of mem:
#$ -l h_vmem=2g

DATA_DIR=$HOME/data/Dmel/H1/
TASKLIST=$DATA_DIR/tasklist.txt

TYPE=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $3}'`
BARCODE=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $4}'`
FILENAME=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $5}'`

# central 95% of combined distribution
MIN_FRAG_LEN=101
MAX_FRAG_LEN=191

INFILE=$DATA_DIR/$FILENAME

echo "$BARCODE / $TYPE" >&2

TRACK_NAME=H1/h1_midpoints_${MIN_FRAG_LEN}_to_${MAX_FRAG_LEN}_$TYPE



CMD="python $HOME/proj/genome/python/script/db/load_solid_mnase_mids.py --assembly dm3 --min_frag_size $MIN_FRAG_LEN --max_frag_size $MAX_FRAG_LEN $TRACK_NAME $INFILE"

echo $CMD >&2
$CMD
