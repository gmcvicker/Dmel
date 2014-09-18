#!/bin/sh
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 1-8
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 4GB of mem:
#$ -l h_vmem=4g

TASKFILE=$HOME/data/Dmel/mnase_midpoint_tracks.txt

LINE=`head -$SGE_TASK_ID $TASKFILE | tail -1`

NAME=`echo $LINE | awk '{print $1}'`
TRACK=`echo $LINE | awk '{print $2}'`

echo "$NAME, hostname=$HOSTNAME" >&2

SCRIPT=$HOME/proj/Dmel/python/mnase_region_counts.py

OUT_DIR=$HOME/data/Dmel/MNase/mnase_region_counts/
mkdir -p $OUT_DIR
OUT_FILE=$OUT_DIR/$NAME.region_counts.txt.gz

python $SCRIPT $TRACK | gzip > $OUT_FILE

if [ "$?" -ne "0" ]; then
    echo "script failed" 1>&2
    exit 1
fi

echo "done" >&2

