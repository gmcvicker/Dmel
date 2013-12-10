#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 2-5
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 8GB of mem:
#$ -l h_vmem=8g

DATA_DIR=$HOME/data/Dmel/HMGD/
TASKLIST=$DATA_DIR/tasklist.txt

REPNAME=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $1}'`
BARCODE=`head -$SGE_TASK_ID $TASKLIST | tail -1 | awk '{print $2}'`

INFILE=$DATA_DIR/bwa_mapped_nb/$BARCODE.BWA.sort.bam

FWD_TRACK=HMGD/hmgd_5prime_ends_rep${REPNAME}_fwd
REV_TRACK=HMGD/hmgd_5prime_ends_rep${REPNAME}_rev
TRACK=HMGD/hmgd_midpoints_rep${REPNAME}

# python $HOME/proj/genome/python/script/db/load_bam_5prime_ends.py --assembly dm3 $FWD_TRACK $REV_TRACK $INFILE

# if [ "$?" -ne "0" ]; then
#     echo "importing BAM file failed" 1>&2
#     exit 1
# fi


# combine fwd/rev tracks, after calculating appropriate offsets for
# each strand
echo "combining ChIP-seq strands" >&2 


python $HOME/proj/script/python/combine_chipseq_strands.py \
    --assembly dm3 $FWD_TRACK $REV_TRACK $TRACK

if [ "$?" -ne "0" ]; then
    echo "combining chipseq strands failed" 1>&2
    exit 1
fi
