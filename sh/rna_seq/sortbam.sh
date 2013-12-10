#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 5-6
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# need 8GB of mem:
#$ -l h_vmem=8g

DIR=/mnt/lustre/home/gmcvicker/data/Dmel/mod_encode/rna_seq

TASK_NAME=`head -$SGE_TASK_ID $DIR/tasklist.txt | tail -1`

BAM_FILE=$DIR/$TASK_NAME.bam

SORTED_BAM_FILE=$DIR/$TASK_NAME.sort.bam

## sort BAM file
echo "sorting BAM" 1>&2
samtools sort -m2000000000 $BAM_FILE $DIR/$TASK_NAME.sort
if [ "$?" -ne "0" ]; then
    echo "samtools sort failed" 1>&2
    exit 1
fi


# index BAM file
echo "indexing BAM file $SORTED_BAM_FILE" 1>&2
samtools index $SORTED_BAM_FILE
if [ "$?" -ne "0" ]; then
    echo "samtools index failed";
    exit 1
fi