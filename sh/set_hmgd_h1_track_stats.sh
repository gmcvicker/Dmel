#!/bin/bash

SCRIPT=$HOME/proj/genome/python/script/db/set_track_stats.py
ASSEMBLY=dm3

for rep in 1 2 3 4 5
do
    echo $rep 1>&2
    python $SCRIPT --assembly $ASSEMBLY HMGD/hmgd_midpoints_rep$rep 
done


for rep in 1 2 3 4
do
    echo $rep 1>&2
    python $SCRIPT --assembly $ASSEMBLY H1/h1_midpoints_101_to_191_rep$rep 
done

