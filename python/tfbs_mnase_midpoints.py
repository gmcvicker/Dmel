import sys
import numpy as np
import re
import os

import genome.db
import genome.coord


HOME = os.environ['HOME']

TFBS_PATH = HOME + "/data/Dmel/redfly/redfly_tfbs_coords.txt"

LINEAGES = ['eye', 'haltere', 'leg', 'antenna']

MNASE_TRACKS = {'eye' : 'mnase/mnase_midpoints_eye_122_to_159',
                'haltere' : 'mnase/mnase_midpoints_haltere_122_to_159',
                'leg' : 'mnase/mnase_midpoints_leg_122_to_159',
                'antenna' : 'mnase/mnase_midpoints_antenna_122_to_159'}

FLANK = 50


class TFBSStat(object):

    def __init__(self):
        self.mnase_count = {}

        for lineage in LINEAGES:
            self.mnase_count[lineage] = 0
        
        self.n_site = 0
        


def get_tfbs(chrom_dict):
    gene_dict = {}

    sys.stderr.write("reading TFBSs from %s\n" % TFBS_PATH)

    f = open(TFBS_PATH, "r")

    tfbs_list = []
    
    for line in f:
        words = line.split()
        chrom = chrom_dict[words[0]]
        start = int(words[1])
        end = int(words[2])
        
        tfbs = genome.coord.Coord(chrom, start, end)

        # first part of TF name is name of factor, second part
        # is name of gene that it is binding near
        # separator is either '-' or '_'
        
        tf_words = words[3].split("_")
        if len(tf_words) < 2:
            raise ValueError("expected '_' separator in TF name '%s'" %
                             words[3])
        tfbs.name = tf_words[0]
        tfbs_list.append(tfbs)

    genome.coord.sort_coords(tfbs_list)

    return tfbs_list



def main():
    gdb = genome.db.GenomeDB(assembly="dm3")
    
    chrom_dict = gdb.get_chromosome_dict()

    tfbs_list = get_tfbs(chrom_dict)

    n_sites = FLANK + FLANK + 1

    mnase_tracks = {}
    for lineage in LINEAGES:
        mnase_tracks[lineage] = gdb.open_track(MNASE_TRACKS[lineage])

    tfbs_stat_dict = {}

    sys.stdout.write("TF.NAME\tCHR\tSTART\tEND\t")
    sys.stdout.write("\t".join([l + ".MNASE.COUNT" for l in LINEAGES]) + "\n")

    # write counts of MNase at each lineage for each TF
    for tfbs in tfbs_list:
        tfbs_mid = tfbs.start + (tfbs.end - tfbs.start + 1)/2
        start = tfbs_mid - FLANK
        end = tfbs_mid + FLANK

        if start < 1:
            continue
        if end > tfbs.chrom.length:
            continue

        sys.stdout.write("%s\t%s\t%d\t%d" %
                         (tfbs.name, tfbs.chrom.name, tfbs.start, tfbs.end))
        
        # get MNase counts for each lineage
        for lineage in LINEAGES:
            track = mnase_tracks[lineage]
            vals = track.get_nparray(tfbs.chrom.name, start, end)
            
            sys.stdout.write("\t%d" % np.sum(vals))


        sys.stdout.write("\n")
    
    for track in mnase_tracks.values():
        track.close()


    


    
main()    
