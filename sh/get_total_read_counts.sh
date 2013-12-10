#!/bin/bash


SCRIPT=$HOME/proj/genome/python/script/db/get_track_stats.py


COUNT=`python $SCRIPT --assembly dm3 mnase/mnase_midpoints_eye_122_to_159 | awk '{print $6}' | awk -F= '{print $2}'`
echo "eye $COUNT" > $HOME/data/Dmel/MNase/mapped_read_counts.txt


COUNT=`python $SCRIPT --assembly dm3 mnase/mnase_midpoints_haltere_122_to_159 | awk '{print $6}' | awk -F= '{print $2}'`
echo "haltere $COUNT" >> $HOME/data/Dmel/MNase/mapped_read_counts.txt

COUNT=`python $SCRIPT --assembly dm3 mnase/mnase_midpoints_leg_122_to_159 | awk '{print $6}' | awk -F= '{print $2}'`
echo "leg $COUNT" >> $HOME/data/Dmel/MNase/mapped_read_counts.txt

COUNT=`python $SCRIPT --assembly dm3 mnase/mnase_midpoints_antenna_122_to_159 | awk '{print $6}' | awk -F= '{print $2}'`
echo "antenna $COUNT" >> $HOME/data/Dmel/MNase/mapped_read_counts.txt

COUNT=`python $SCRIPT --assembly dm3 H1/h1_midpoints_101_to_191_combined | awk '{print $6}' | awk -F= '{print $2}'`
echo "h1 $COUNT" >> $HOME/data/Dmel/MNase/mapped_read_counts.txt

COUNT=`python $SCRIPT --assembly dm3 HMGD/hmgd_midpoints_combined | awk '{print $6}' | awk -F= '{print $2}'`
echo "hmgd $COUNT" >> $HOME/data/Dmel/MNase/mapped_read_counts.txt

COUNT=`python $SCRIPT --assembly dm3 mnase/mnase_midpoints_S2_101_to_191 | awk '{print $6}' | awk -F= '{print $2}'`
echo "S2 $COUNT" >> $HOME/data/Dmel/MNase/mapped_read_counts.txt


COUNT=`python $SCRIPT --assembly dm3 dnase/dnase_S2 | awk '{print $6}' | awk -F= '{print $2}'`
echo "dnase $COUNT" >> $HOME/data/Dmel/MNase/mapped_read_counts.txt


