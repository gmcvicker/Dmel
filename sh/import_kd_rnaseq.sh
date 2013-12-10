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
# jobs need 4GB of mem:
#$ -l h_vmem=4g

DATA_DIR=/home/gmcvicker/data/Dmel/KD_RNA_seq/
TASKLIST=$DATA_DIR/tasklist.txt

TYPE=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $2}'`

BAMFILE_WITH_DUPS=$DATA_DIR/bwa_mapped_nb/${TYPE}_1.BWA.sort.bam
BAMFILE=$DATA_DIR/bwa_mapped_nb/${TYPE}_1.BWA.sort.rmdup.bam

TRACK_NAME=rna_seq/S2_KD/${TYPE}
FWD_TRACK_NAME=${TRACK_NAME}_fwd
REV_TRACK_NAME=${TRACK_NAME}_rev

echo "removing duplicate reads" >&2
CMD="samtools rmdup -s $BAMFILE_WITH_DUPS $BAMFILE"
echo "CMD=$CMD" >&2
$CMD
if [ "$?" -ne "0" ]; then
    echo "samtools rmdup failed" 1>&2
    exit 1
fi


echo "indexing bam file" >&2
samtools index $BAMFILE

if [ "$?" -ne "0" ]; then
    echo "samtools index failed" 1>&2
    exit 1
fi


echo "$TYPE" >&2
echo "importing read counts from BAM" >&2
IMPORT_SCRIPT=/home/gmcvicker/proj/genome/python/script/db/load_bam_5prime_ends.py
python $IMPORT_SCRIPT --assembly dm3 $FWD_TRACK_NAME $REV_TRACK_NAME $BAMFILE

if [ "$?" -ne "0" ]; then
    echo "load_bam_5prime_ends failed" 1>&2
    exit 1
fi


echo "combining fwd/rev tracks" >&2
COMBINE_SCRIPT=/home/gmcvicker/proj/genome/python/script/db/combine_tracks.py
python $COMBINE_SCRIPT --assembly dm3 --dtype uint8 $TRACK_NAME $FWD_TRACK_NAME $REV_TRACK_NAME

if [ "$?" -ne "0" ]; then
    echo "combine_tracks failed" 1>&2
    exit 1
fi

echo "setting track stats" >&2
STAT_SCRIPT=/home/gmcvicker/proj/genome/python/script/db/set_track_stats.py
python $STAT_SCRIPT --assembly dm3  $TRACK_NAME

if [ "$?" -ne "0" ]; then
    echo "set_track_stats failed" 1>&2
    exit 1
fi

python $STAT_SCRIPT --assembly dm3  $REV_TRACK_NAME

if [ "$?" -ne "0" ]; then
    echo "set_track_stats (rev) failed" 1>&2
    exit 1
fi

python $STAT_SCRIPT --assembly dm3  $FWD_TRACK_NAME

if [ "$?" -ne "0" ]; then
    echo "set_track_stats (fwd) failed" 1>&2
    exit 1
fi
