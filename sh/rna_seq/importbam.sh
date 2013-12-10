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

DIR=$HOME/data/Dmel/mod_encode/rna_seq

SCRIPT_DIR=$HOME/proj/genome/python/script/db
COMBINE_SCRIPT=$SCRIPT_DIR/combine_tracks.py
STAT_SCRIPT=$SCRIPT_DIR/set_track_stats.py

ASSEMBLY=dm3

TASK_NAME=`head -$SGE_TASK_ID $DIR/tasklist.txt | tail -1`

BAM_FILE=$DIR/$TASK_NAME.bam

SORTED_BAM_FILE=$DIR/$TASK_NAME.sort.bam

TRACK_NAME=`echo $TASK_NAME | sed s/-/_/g`


# for TYPE in 5prime_ends read_depth
for TYPE in read_depth
do
    COMBINED_TRACK=rna_seq/${TRACK_NAME}_${TYPE}
    IMPORT_SCRIPT=$SCRIPT_DIR/load_bam_$TYPE.py
	

    if [ "$TYPE" == "5prime_ends" ]; then
    
	FWD_TRACK=rna_seq/${TRACK_NAME}_${TYPE}_fwd
	REV_TRACK=rna_seq/${TRACK_NAME}_${TYPE}_rev

	python $IMPORT_SCRIPT --assembly $ASSEMBLY \
	    $FWD_TRACK $REV_TRACK $SORTED_BAM_FILE

	if [ "$?" -ne "0" ]; then
	    echo "import failed" 1>&2
	    exit 1
	fi

	python $COMBINE_SCRIPT --assembly $ASSEMBLY \
	    $COMBINED_TRACK $FWD_TRACK $REV_TRACK

	if [ "$?" -ne "0" ]; then
	    echo "combine tracks failed" 1>&2
	    exit 1
	fi

	python $STAT_SCRIPT --assembly $ASSEMBLY $FWD_TRACK
	
	if [ "$?" -ne "0" ]; then
	    echo "setting stats for fwd track $FWD_TRACK failed" 1>&2
	    exit 1
	fi

	python $STAT_SCRIPT --assembly $ASSEMBLY $REV_TRACK
	
	if [ "$?" -ne "0" ]; then
	    echo "setting stats for rev track $REV_TRACK failed" 1>&2
	    exit 1
	fi
    else
	python $IMPORT_SCRIPT --assembly $ASSEMBLY \
	    $COMBINED_TRACK $SORTED_BAM_FILE
    fi

    python $STAT_SCRIPT --assembly $ASSEMBLY $COMBINED_TRACK

    if [ "$?" -ne "0" ]; then
	echo "setting stats for rev track $COMBINED_TRACK failed" 1>&2
	exit 1
    fi
done

