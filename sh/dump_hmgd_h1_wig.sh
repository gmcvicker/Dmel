
SCRIPT=$HOME/proj/script/python/dump_wig.py

OUTPUT_DIR=$HOME/data/Dmel/HMGD_H1_export
GEO_DIR=$HOME/data/Dmel/HMGD_H1_GEO_SUBMISSION

TRACKS="HMGD/hmgd_midpoints_rep1 HMGD/hmgd_midpoints_rep2 HMGD/hmgd_midpoints_rep3 HMGD/hmgd_midpoints_rep5 H1/h1_midpoints_101_to_191_rep2 H1/h1_midpoints_101_to_191_rep4 mnase/mnase_midpoints_S2_101_to_191"

# TRACKS="HMGD/hmgd_midpoints_rep5 mnase/mnase_midpoints_S2_101_to_191"

for TRACK in $TRACKS
do
    # echo $TRACK
    PREFIX=`echo $TRACK | sed 's/\//_/g'`
    echo $PREFIX
    mkdir -p $OUTPUT_DIR/$PREFIX
    # python $SCRIPT --assembly dm3 --combine_files $TRACK $OUTPUT_DIR/$PREFIX

    cp $OUTPUT_DIR/$PREFIX/combined.wig.gz $GEO_DIR/$PREFIX.wig.gz
done

