

SCRIPT=/home/gmcvicker/proj/script/python/Dmel/find_tss_motif.py

OUTDIR=/home/gmcvicker/data/Dmel/tss_motif/

mkdir -p $OUTDIR


# echo "CCAAT" >&2
# python $SCRIPT --max_mismatch 1 -80 -65 CCAAT > $OUTDIR/ccaat.txt

# echo "TATA box" >&2
# python $SCRIPT --max_mismatch 1 -40 -20 TATAWAWR > $OUTDIR/tatabox.txt

echo "TATA box" >&2
python $SCRIPT --max_mismatch 1 -40 -20 RWAWATAT > $OUTDIR/tatabox.rev.txt

# echo "TATA box" >&2
# python $SCRIPT --max_mismatch 1 -400 -380 TATAWAWR > $OUTDIR/tatabox.move.txt

# echo "BRE_U" >&2
# python $SCRIPT --max_mismatch 1 -50 -20 SSRCGCC > $OUTDIR/bre.u.txt

# echo "BRE_D" >&2
# python $SCRIPT --max_mismatch 1 -40 -20 RTDKKKK > $OUTDIR/bre.d.txt

# echo "Inr" >&2
# python $SCRIPT --max_mismatch 1 -10 10 TCAKTY > $OUTDIR/inr.txt

# echo "DPE" >&2
# python $SCRIPT --max_mismatch 1 25 40 RGWYV > $OUTDIR/dpe.txt

# echo "MTE" >&2
# python $SCRIPT --max_mismatch 1 15 30 CSARCSSAAC > $OUTDIR/mte.txt

