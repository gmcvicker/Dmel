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
# jobs need 4GB of mem:
#$ -l h_vmem=4g

TASKFILE=$HOME/data/Dmel/lineage.txt

LINE=`head -$SGE_TASK_ID $TASKFILE | tail -1`

CELL_LINE=`echo $LINE | awk '{print $2}'`
LINEAGE_NAME=`echo $LINE | awk '{print $3}'`
BARCODE=`echo $LINE | awk '{print $4}'`

KMER_SIZE=7

echo "$CELL_LINE($LINEAGE_NAME), hostname=$HOSTNAME" >&2


TRACK=mnase/mnase_midpoints_${LINEAGE_NAME}_122_to_159
SCRIPT=$HOME/proj/script/python/Dmel/mnase_kmer_counts.py
OUT_FILE=$HOME/data/Dmel/MNase/kmer_counts/$LINEAGE_NAME.${KMER_SIZE}mer.txt

python $SCRIPT $KMER_SIZE $TRACK > $OUT_FILE

if [ "$?" -ne "0" ]; then
    echo "script failed" 1>&2
    exit 1
fi

echo "done" >&2

