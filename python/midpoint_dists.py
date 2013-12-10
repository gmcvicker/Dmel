
import sys

import genome.db

import util.sample
import numpy as np


LINEAGES = ['eye', 'haltere', 'leg', 'antenna']

MNASE_TRACKS = {'eye' : 'mnase/mnase_midpoints_eye_122_to_159',
                'haltere' : 'mnase/mnase_midpoints_haltere_122_to_159',
                'leg' : 'mnase/mnase_midpoints_leg_122_to_159',
                'antenna' : 'mnase/mnase_midpoints_antenna_122_to_159'}


MAX_DIST = 10000

gdb = genome.db.GenomeDB(assembly='dm3')


dist_counts = {}
tracks = {}
for lineage in LINEAGES:
    dist_counts[lineage] = np.zeros(MAX_DIST+1, dtype=np.int32)
    tracks[lineage] = gdb.open_track(MNASE_TRACKS[lineage])

for chrom in gdb.get_chromosomes():
    sys.stderr.write("%s\n" % chrom.name)
    
    for lineage in LINEAGES:
        sys.stderr.write("  %s\n" % lineage)
        
        vals = tracks[lineage].get_nparray(chrom)

        indices = np.where(vals > 0)[0]

        samp_size = 1e4
        sys.stderr.write("  sampling %d sites\n" % samp_size)
        
        samp = util.sample.weighted_sample(indices, vals[indices], samp_size)
        
        for idx in samp:
            if idx + MAX_DIST >= chrom.length:
                continue
            if idx - MAX_DIST <= 0:
                continue

            count_at_site = vals[idx]
            right_counts = vals[(idx+1):(idx+MAX_DIST+1)]
            left_counts = vals[(idx-1):(idx-MAX_DIST-1):-1]
            dist_counts[lineage][1:] += right_counts * count_at_site
            dist_counts[lineage][1:] += left_counts * count_at_site
            dist_counts[lineage][0] += count_at_site
            

# write header
sys.stdout.write("DIST " + " ".join(LINEAGES) + "\n")

# write out counts for each distance
for i in range(MAX_DIST+1):
    sys.stdout.write("%d" % i)
    for lineage in LINEAGES:
        sys.stdout.write(" %d" % dist_counts[lineage][i])
    sys.stdout.write("\n")


for track in tracks.values():
    track.close()
