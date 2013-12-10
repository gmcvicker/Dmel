#!/bin/bash

SCRIPT=$HOME/proj/genome/python/script/db/set_track_stats.py
ASSEMBLY=dm3

for LINEAGE in antenna eye haltere leg
do
    echo $LINEAGE 1>&2
    python $SCRIPT --assembly $ASSEMBLY mnase/mnase_midpoints_${LINEAGE}_122_to_159
done

