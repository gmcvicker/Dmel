#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 2GB of mem:
#$ -l h_vmem=2g

# this is set to match imported H1, since we are using MNase as background
# for H1
MIN_FRAG_LEN=101
MAX_FRAG_LEN=191

INFILE=$HOME/data/Dmel/MNase-S2/WIDOM_061311RLM_PE_S2_Viv_AACATG_paired_results.txt.gz

echo "$BARCODE / $TYPE" >&2

TRACK_NAME=mnase/mnase_midpoints_S2_${MIN_FRAG_LEN}_to_${MAX_FRAG_LEN}

CMD="python $HOME/proj/genome/python/script/db/load_solid_mnase_mids.py --assembly dm3 --min_frag_size $MIN_FRAG_LEN --max_frag_size $MAX_FRAG_LEN $TRACK_NAME $INFILE"

echo $CMD >&2
$CMD
