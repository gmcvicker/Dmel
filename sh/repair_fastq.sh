#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 1-5
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 8GB of mem:
#$ -l h_vmem=8g


DATA_DIR=$HOME/data/Dmel/HMGD/
TASKFILE=$DATA_DIR/tasklist.txt

LINE=`head -$SGE_TASK_ID $TASKFILE | tail -1`
BARCODE=`echo $LINE | awk '{print $2}'`

# input file
FASTQ_FILE=$DATA_DIR/fastq/YFM_${BARCODE}_061212.fastq.gz

# output file
OUTPUT_FILE=$DATA_DIR/fastq/YFM_${BARCODE}_061212.repaired.fastq.gz

# repair file
echo "repairing file $FASTQ_FILE" 1>&2
$HOME/proj/check_fastq/check_fastq $FASTQ_FILE $OUTPUT_FILE 
if [ "$?" -ne "0" ]; then
    echo "check_fastq failed";
    exit 1
fi


echo "done" 1>&2

