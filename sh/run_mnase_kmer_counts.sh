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

KMER_SIZE=7

echo "$NAME, hostname=$HOSTNAME" >&2

SCRIPT=$HOME/proj/Dmel/python/mnase_kmer_counts.py
OUT_FILE=$HOME/data/Dmel/MNase/kmer_counts/$NAME.${KMER_SIZE}mer.txt

python $SCRIPT $KMER_SIZE $TRACK > $OUT_FILE

if [ "$?" -ne "0" ]; then
    echo "script failed" 1>&2
    exit 1
fi

echo "done" >&2

