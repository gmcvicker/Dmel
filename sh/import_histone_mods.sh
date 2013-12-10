#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
#$ -t 1-102
#
# jobs need 2GB of mem, quite a bit of IO:
#$ -l h_vmem=2g,bigio=1

TASKLIST=$HOME/data/Dmel/mod_encode/histone_modification/histone_mod_tasklist.txt

MARK=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $2}'`
CELL_LINE=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $3}'`
EXPERIMENT_ID=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $5}'`
FILENAME=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $6}'`

#
# Split files by chromosome
# note files are actually similar to bedgraph format
# even though they end with 'wig'
#
TMP_FILENAME=${MARK}_${CELL_LINE}_${EXPERIMENT_ID}.bedgraph.gz
ln -s $FILENAME /tmp/$TMP_FILENAME
CMD="$HOME/proj/genome/c/program/split_bed_chrs /tmp/$TMP_FILENAME"
echo $CMD >&2
$CMD

if [ "$?" -ne "0" ]; then
    echo "splitting chromosomes failed" 1>&2
    exit 1
fi

#
# create new track
#
echo "$MARK $CELL_LINE $EXPERIMENT_ID" >&2

TRACK_NAME=histone_mods/${MARK}_${CELL_LINE}_${EXPERIMENT_ID}
INPUT_FILENAMES=`ls /tmp/chr*_$TMP_FILENAME | xargs`
CMD="python $HOME/proj/genome/python/script/db/create_track.py --format bedgraph --dtype float32 --assembly dm3 $TRACK_NAME $INPUT_FILENAMES"
echo $CMD >&2
$CMD

if [ "$?" -ne "0" ]; then
    echo "import failed" 1>&2
    exit 1
fi

#
# set track stats
#
CMD="python $HOME/proj/genome/python/script/db/set_track_stats.py --assembly dm3 $TRACK_NAME"
echo $CMD >&2
$CMD

if [ "$?" -ne "0" ]; then
    echo "set_track_stats failed" 1>&2
    exit 1
fi


#
# cleanup temporary files
#
rm /tmp/chr*_$TMP_FILENAME
rm /tmp/$TMP_FILENAME

echo "done" >&2

