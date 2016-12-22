
python ../python/tss_mnase_midpoints.py | gzip > ~/data/Dmel/tss_mnase_midpoints.all.txt.gz
python ../python/tss_mnase_midpoints.py --low_expr | gzip > ~/data/Dmel/tss_mnase_midpoints.low_expr.txt.gz
python ../python/tss_mnase_midpoints.py --mid_expr | gzip > ~/data/Dmel/tss_mnase_midpoints.mid_expr.txt.gz
python ../python/tss_mnase_midpoints.py --high_expr | gzip > ~/data/Dmel/tss_mnase_midpoints.high_expr.txt.gz
