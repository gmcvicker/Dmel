#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 1-3
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 8GB of mem:
#$ -l h_vmem=8g


INDEX_FILE=$HOME/bwa/bwa_indexes_nb/dm3.fa
DATA_DIR=$HOME/data/Dmel/DNase
TASKFILE=$DATA_DIR/tasklist.txt

OUT_DIR=$DATA_DIR/bwa_mapped_nb

SAMTOOLS=/usr/local/bin/samtools

# input file
FILENAME=`head -$SGE_TASK_ID $TASKFILE | tail -1`
FASTQ_FILE=$DATA_DIR/$FILENAME
PREFIX=`echo $FILENAME | sed "s/.fastq.gz$//"`

echo "hostname: $HOSTNAME filename: $FASTQ_FILE" >&2

# output files 
SAI_FILE=$OUT_DIR/$PREFIX.BWA.sai
SAM_FILE=$OUT_DIR/$PREFIX.BWA.sam.gz
BAM_FILE=$OUT_DIR/$PREFIX.BWA.bam
SORTED_BAM_PREFIX=$OUT_DIR/$PREFIX.BWA.sort
SORTED_BAM_FILE=$SORTED_BAM_PREFIX.bam

mkdir -p $OUT_DIR

# Align reads, output suffix array (SA) coordinates
echo "aligning reads" 1>&2
echo "command: bwa aln $INDEX_FILE $FASTQ_FILE > $SAI_FILE" >&2
bwa aln $INDEX_FILE $FASTQ_FILE > $SAI_FILE
if [ "$?" -ne "0" ]; then
    echo "bwa aln failed" 1>&2
    exit 1
fi

# convert SA coordinates to SAM alignment format
echo "generating SAM file" 1>&2
bwa samse $INDEX_FILE $SAI_FILE $FASTQ_FILE | gzip > $SAM_FILE
if [ "$?" -ne "0" ]; then
    echo "bwa sampe failed" 1>&2
    exit 1
fi

# # rm $SAI_FILE


# convert SAM to binary BAM format
echo "creating BAM file" 1>&2
$SAMTOOLS view -S -b $SAM_FILE > $BAM_FILE
if [ "$?" -ne "0" ]; then
    echo "samtools view failed" 1>&2
    exit 1
fi

## OLD version of samtools:
## samtools import $INDEX_FILE $SAM_FILE $BAM_FILE

# rm $SAM_FILE

## sort BAM file
echo "sorting BAM" 1>&2
$SAMTOOLS sort -m2000000000 $BAM_FILE $SORTED_BAM_PREFIX
if [ "$?" -ne "0" ]; then
    echo "samtools sort failed" 1>&2
    exit 1
fi

# rm $BAM_FILE


# index BAM file
echo "indexing BAM file $SORTED_BAM_FILE" 1>&2
$SAMTOOLS index $SORTED_BAM_FILE
if [ "$?" -ne "0" ]; then
    echo "samtools index failed";
    exit 1
fi


echo "done" 1>&2

