#!/bin/bash

SCRIPT=$HOME/proj/Dmel/python/kmers/region_kmer_counts.py

DATA_DIR=$HOME/data/Dmel/MNase/mnase_region_counts/signif_diff/

for LINEAGE in S2_in_vitro_combined antenna eye haltere leg
do
    echo $LINEAGE 1>&2
    PREFIX=$DATA_DIR/${LINEAGE}_fdr0.01
    python $SCRIPT ${PREFIX}_higher.txt.gz > ${PREFIX}_higher_kmers.txt
    python $SCRIPT ${PREFIX}_lower.txt.gz > ${PREFIX}_lower_kmers.txt
    
done

