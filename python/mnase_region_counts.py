import sys
import numpy as np
import re
import os

import genome.db
import genome.coord


HOME = os.environ['HOME']

LINEAGES = ['eye', 'haltere', 'leg', 'antenna']

MNASE_TRACKS = {'eye' : 'mnase/mnase_midpoints_eye_122_to_159',
                'haltere' : 'mnase/mnase_midpoints_haltere_122_to_159',
                'leg' : 'mnase/mnase_midpoints_leg_122_to_159',
                'antenna' : 'mnase/mnase_midpoints_antenna_122_to_159'}

REGION_SIZE = 200


            

def get_total_counts(gdb, mnase_track, other_mnase_tracks):
    stat = gdb.get_track_stat(mnase_track)

    total = stat.sum
    
    other_total = 0
    for other_track in other_mnase_tracks:
        other_stat = gdb.get_track_stat(other_track)
        other_total += other_stat.sum

    return (total, other_total)
    
    
    

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <lineage>\n" % sys.argv[0])
        exit(1)
    
    gdb = genome.db.GenomeDB(assembly='dm3')

    # open mnase count tracks for the 'main' lineage and the 'other' lineages
    # that we are comparing to
    lineage = sys.argv[1]
    if lineage not in LINEAGES:
        raise ValueError("lineage must be one of %s" % ",".join(LINEAGES))

    mnase_track = gdb.open_track(MNASE_TRACKS[lineage])
    other_mnase_tracks = []
    for lin in LINEAGES:
        if lin != lineage:
            other_mnase_tracks.append(gdb.open_track(MNASE_TRACKS[lin]))


    (total_counts, other_total_counts) = get_total_counts(gdb, mnase_track,
                                                          other_mnase_tracks)
                                                          

    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom.name)

        start = 1
        end = start + REGION_SIZE - 1

        # divide genome into non-overlapping regions
        while start < chrom.length:            
            # get count of midpoints for main lineage and combined
            # count for other lineages for this region
            
            main_counts = np.sum(mnase_track.get_nparray(chrom, start, end))

            other_counts = 0
            for track in other_mnase_tracks:
                other_counts += np.sum(track.get_nparray(chrom, start, end))
                
            sys.stdout.write("%s %d %d %d %d %d %d\n" % 
                             (chrom.name, start, end,
                              main_counts, total_counts,
                              other_counts, other_total_counts))
                                    
            start += REGION_SIZE
            end = start + REGION_SIZE - 1

            if end > chrom.length:
                end = chrom.length

            

            
main()
            


        
    
