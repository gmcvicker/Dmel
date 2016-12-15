
for CELLTYPE in antenna eye haltere leg S2; do
    echo $CELLTYPE >& 2
    python ../python/mnase_midpoint_profile.py mnase/mnase_midpoints_${CELLTYPE}_122_to_159 ~/data/Dmel/${CELLTYPE}_profile.txt
done


