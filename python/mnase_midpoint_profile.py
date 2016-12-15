import sys
import os.path

import argparse
import numpy as np

import util.sample
import genome.db
import profile


NUCSOME_SIZE = 147
NUCSOME_MID = 74
FLANKING = 5000

MAX_WEIGHT = 20

    


def get_positions(gdb, options, chrom, weights, tracks, profile_len):
    # set beginning and end of chrs to 0s so they are not used
    weights[:profile_len] = 0
    weights[-profile_len:] = 0

    f = (weights > 0.0)

    return np.where(f)[0] + 1
    
    

def parse_arguments():
    parser = argparse.ArgumentParser(description="create aggregate profiles "
                                     "of MNase, sequence composition "
                                     "etc. around putative nucleosome dyad "
                                     "positions")

    parser.add_argument("track_name", help="name of genome DB track to obtain "
                       "midpoints from")

    parser.add_argument('output_file', help='file to write output to')


    parser.add_argument("--chrom", help="run on specified chromosome",
                        default=None)
    
    args = parser.parse_args()

    if os.path.exists(args.output_file):
        raise IOError("output file '%s' already exists\n" % args.output_file)
    else:
        sys.stderr.write("writing output to file '%s'\n" % args.output_file)

    return args
    


def main():
    args = parse_arguments()
    
    output_file = open(args.output_file, "w")

    gdb = genome.db.GenomeDB(assembly="dm3")
    
    chromosomes = gdb.get_chromosomes()

    if args.chrom:
        chrom_dict = gdb.get_chromosome_dict()
        if args.chrom in chrom_dict:
            chromosomes = [chrom_dict[args.chrom]]
        else:
            ValueError("Unknown chromosome '%s'. Chromosomes: %s\n" %
                       (args.chrom, ",".join(chromosomes)))
    
    
    tracks = profile.DataTracks(gdb, args.track_name)
    weight_track = tracks.mnase
        
    profile_len = FLANKING + FLANKING + NUCSOME_SIZE
    profile_mid = FLANKING + NUCSOME_MID
    pfl = profile.Profile(profile_len, profile_mid, output_file)

    count = 0
    
    for chrom in chromosomes:
        sys.stderr.write("%s\n" % chrom.name)
        sys.stderr.write("reading data\n")

        # threshold weights
        weights = weight_track.get_nparray(chrom)
        weights[weights > MAX_WEIGHT] = MAX_WEIGHT
                        
        pos_list = get_positions(gdb, args, chrom, weights,
                                 tracks, profile_len)
        ttl_weight = np.sum(weights)

        sys.stderr.write("  %d positions, total weight=%g\n" %
                         (pos_list.size, ttl_weight))
        
        sys.stderr.write("processing data\n")
        count = 0
        for pos in pos_list:
            weight = weights[pos-1]

            count += 1
            if count > 1000:
                sys.stderr.write(".")
                count = 0
            
            start = pos - NUCSOME_MID - FLANKING + 1
            end = pos + NUCSOME_MID + FLANKING - 1
            pfl.add_counts(tracks, chrom, start, end, 1, weight=weight)
            
        sys.stderr.write("\n")

    pfl.write()
    tracks.close()
    weight_track.close()
    pfl.out_f.close()
    

if __name__ == "__main__":
    main()
