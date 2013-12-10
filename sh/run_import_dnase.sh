#!/bin/sh

DATA_DIR=$HOME/data/Dmel/DNase/bwa_mapped_nb

SCRIPT_DIR=$HOME/proj/genome/python/script/db
LOAD_SCRIPT=$SCRIPT_DIR/load_bam_5prime_ends.py
MERGE_SCRIPT=$SCRIPT_DIR/combine_tracks.py
STAT_SCRIPT=$SCRIPT_DIR/set_track_stats.py


FWD_TRACK=dnase/dnase_S2_fwd
REV_TRACK=dnase/dnase_S2_rev
COMBINED_TRACK=dnase/dnase_S2

echo "loading read counts" >&2
python $LOAD_SCRIPT --assembly dm3 $FWD_TRACK $REV_TRACK $DATA_DIR/S2_DS15591_FC6228W_6.BWA.sort.bam $DATA_DIR/S2_DS15592_FC6228W_7.BWA.sort.bam $DATA_DIR/S2_FC42YCR-6.BWA.sort.bam

echo "merging tracks" >&2
python $MERGE_SCRIPT --assembly dm3  $COMBINED_TRACK $FWD_TRACK $REV_TRACK

echo "setting track stats" >&2
python $STAT_SCRIPT --assembly dm3 $FWD_TRACK 
python $STAT_SCRIPT --assembly dm3 $REV_TRACK
python $STAT_SCRIPT --assembly dm3 $COMBINED_TRACK

