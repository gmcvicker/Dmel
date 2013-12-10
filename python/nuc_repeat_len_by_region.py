import sys


import genome.db
import genome.coord

import numpy as np
import gzip


MNASE_TRACK = "mnase/mnase_midpoints_S2_101_to_191"

MAX_DIST = 3000



def read_regions(filename, chrom_dict):    
    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)

    header = f.readline()
    
    coords = dict([(x, []) for x in chrom_dict.keys()])
    
    for line in f:
        words = line.rstrip().split()
        chrom = chrom_dict[words[0]]
        start = int(words[1])
        end = int(words[2])

        coord = genome.coord.Coord(chrom, start, end)
        
        coords[chrom.name].append(coord)

    f.close()
        
    return coords
    

    
def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: %s <regions.bed.gz> <min_count> > "
                         "output.txt\n" % sys.argv[0])
        exit(2)

    input_filename = sys.argv[1]
    min_count = int(sys.argv[2])
        
    gdb = genome.db.GenomeDB(assembly="dm3")

    mnase_track = gdb.open_track(MNASE_TRACK)
    
    sys.stdout.write("DIST N.SITES COUNT\n")

    chrom_dict = gdb.get_chromosome_dict()
    coord_dict = read_regions(input_filename, chrom_dict)

    counts_by_dist = np.zeros(MAX_DIST+1, dtype=np.int64)
    sites_by_dist = np.zeros(MAX_DIST+1, dtype=np.int64)
    n_sites = 0
    
    for chrom in gdb.get_chromosomes():

        coords = coord_dict[chrom.name]

        vals = mnase_track.get_nparray(chrom)
        
        sys.stderr.write("%s\n" % chrom.name)

        for coord in coords:

            # get index of where there are >= min_count midpoints 
            # in this region
            idx = np.where(vals[coord.start-1:coord.end] >= min_count)[0]
            idx += coord.start - 1

            if coord.start - MAX_DIST < 1:
                continue
            if coord.end + MAX_DIST > chrom.length:
                continue

            for i in idx:
                # get counts to left of this midpoint 
                # (not including midpoint)
                start = i - MAX_DIST
                end   = i 
                left_counts = vals[start:end][::-1]
                
                # get counts to right of this midpoint
                start = i + 1
                end   = i + MAX_DIST + 1
                right_counts = vals[start:end]

                mid_count = vals[i]
                counts_by_dist[1:] += left_counts * mid_count
                counts_by_dist[1:] += right_counts * mid_count
                counts_by_dist[0] += mid_count

                sites_by_dist[1:] += mid_count * 2
                sites_by_dist[0] += mid_count


    for i in range(MAX_DIST+1):
        sys.stdout.write("%d %d %d\n" % (i, sites_by_dist[i], counts_by_dist[i]))

    mnase_track.close()
    



main()
