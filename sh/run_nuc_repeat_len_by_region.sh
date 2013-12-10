#!/bin/bash

SCRIPT=$HOME/proj/script/python/Dmel/nuc_repeat_len_by_region.py
DIR=$HOME/data/Dmel/nuc_repeat_len_by_region/

for TYPE in hmgd_h1_ratio h1 hmgd
do
    for Q in {1..5}
    do
	echo ${TYPE}_Q$Q
	python $SCRIPT $DIR/regions_${TYPE}_Q$Q.txt 3 > $DIR/nrl_${TYPE}_Q$Q.txt
    done
done

