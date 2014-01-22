#!/bin/sh
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
# jobs need 2GB of mem
#$ -l h_vmem=2g


SCRIPT=$HOME/proj/genome/python/script/db/load_solid_mnase_mids.py

# central 95% of frag sizes from D. melanogaster samples
# MIN_FRAG_SIZE=122
# MAX_FRAG_SIZE=159

# match sizes used for H1 and HMGD1 samples
MIN_FRAG_SIZE=101
MAX_FRAG_SIZE=191

DATA_DIR=$HOME/data/Dmel/MNase
TASKLIST=$DATA_DIR/tasklist.txt

TYPE=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $3}'`
BARCODE=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $4}'`
FILENAME=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $5}'`

INFILE=$DATA_DIR/$BARCODE/$FILENAME

TRACK_NAME=/mnase/mnase_midpoints_${TYPE}_${MIN_FRAG_SIZE}_to_${MAX_FRAG_SIZE}

echo "$BARCODE / $TYPE" >&2
echo "python $SCRIPT --assembly dm3 --min_frag_size $MIN_FRAG_SIZE --max_frag_size $MAX_FRAG_SIZE $TRACK_NAME $INFILE" >&2

python $SCRIPT --assembly dm3 --min_frag_size $MIN_FRAG_SIZE --max_frag_size $MAX_FRAG_SIZE $TRACK_NAME $INFILE

# set track stats
STAT_SCRIPT=$HOME/proj/genome/python/script/db/set_track_stats.py
echo "python $STAT_SCRIPT --assembly dm3 $TRACK_NAME" >&2
python $STAT_SCRIPT --assembly dm3 $TRACK_NAME
