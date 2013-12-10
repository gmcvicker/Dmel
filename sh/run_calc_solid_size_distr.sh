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

INFILE=$DATA_DIR/$FILENAME
OUTFILE=$DATA_DIR/frag_size_distr/$TYPE.txt

echo "$BARCODE / $TYPE" >&2
echo "python $HOME/proj/script/python/calc_solid_frag_size_distribution.py --assembly dm3 $INFILE > $OUTFILE" >&2

python $HOME/proj/script/python/calc_solid_frag_size_distribution.py --assembly dm3 $INFILE > $OUTFILE

