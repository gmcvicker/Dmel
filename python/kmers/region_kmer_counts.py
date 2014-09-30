import sys
import numpy as np
import re
import os

import genome.kmer

import genome.db
import genome.coord
import gzip

KMER_SIZE=6


#
# Counts kmers occurances in a specified set of genomic regions, as well as genome-wide
#



def read_regions(filename, chrom_dict):
    coord_dict = {}
    for chrom_name in chrom_dict.keys():
        coord_dict[chrom_name] = []
        

    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)

    # skip header line if there is one
    first_line = f.readline()
    if not first_line.startswith("CHROM"):
        # this is not a header, rewind file
        f.seek(0)

    for line in f:
        words = line.rstrip().split()
        chrom = chrom_dict[words[0]]
        start = int(words[1])
        end = int(words[2])

        coord = genome.coord.Coord(chrom, start, end)
        coord_dict[chrom.name].append(coord)

    return coord_dict
                           
        
    
    

def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <region_file>\n" % sys.argv[0])
        exit(1)
    
    gdb = genome.db.GenomeDB(assembly='dm3')

    chrom_dict = gdb.get_chromosome_dict()
    coord_dict = read_regions(sys.argv[1], chrom_dict)

    seq_track = gdb.open_track("seq")


    all_kmers = {}
    obs_kmers = {}
    
    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom)
        
        # flag regions from coordinate file        
        sys.stderr.write("  flagging regions\n")
        coords = coord_dict[chrom.name]
        flags = np.zeros(chrom.length, np.int32)
        for coord in coords:
            flags[(coord.start-1):coord.end] = 1

        sys.stderr.write("  getting DNA sequence\n")
        dna = seq_track.get_seq_str(chrom)

        sys.stderr.write("  counting kmers\n")
        # count kmers on this chromosome
        genome.kmer.count_kmers(dna, KMER_SIZE, flags, all_kmers, obs_kmers)
        
    sys.stdout.write("KMER REGION.COUNT ALL.COUNT\n")
    for kmer_str in genome.kmer.list_kmers(KMER_SIZE):
        if kmer_str in obs_kmers:
            obs_count = obs_kmers[kmer_str]
        else:
            obs_count = 0

        if kmer_str in all_kmers:
            all_count = all_kmers[kmer_str]
        else:
            obs_count = 0
            
        sys.stdout.write("%s %d %d\n" % (kmer_str, obs_count, all_count))


    seq_track.close()
        

main()
            


        
    
