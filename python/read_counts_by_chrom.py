import genome.db
import genome.trackstat
import sys


gdb = genome.db.GenomeDB(assembly='dm3')

track_name = sys.argv[1]
track = gdb.open_track(track_name)

sys.stdout.write("CHROM\tCHROM.LEN\tN.MIDPOINTS\tMIDPOINTS.PER.SITE\n")
for chrom in gdb.get_chromosomes():
    stat = genome.trackstat.get_stats(gdb, track, chrom)
    mids_per_site = float(stat.sum)/chrom.length
    sys.stdout.write("%s\t%d\t%d\t%.4f\n" % 
                     (chrom.name, chrom.length, stat.sum, mids_per_site))

track.close()
