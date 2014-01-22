import sys
import numpy as np
import re
import os

import genome.kmer

import genome.db
import genome.coord

MNASE_SMOOTH_WINDOW = 150

def main():
    if len(sys.argv) != 3:
        sys.stderr.write("usage: %s <kmer_size> <mnase_midpoint_track>\n" % sys.argv[0])
        exit(1)
    
    gdb = genome.db.GenomeDB(assembly='dm3')

    seq_track = gdb.open_track("seq")

    kmer_size = int(sys.argv[1])
    mnase_track = gdb.open_track(sys.argv[2])
    
    all_kmers = {}
    obs_kmers = {}

    window = np.ones(MNASE_SMOOTH_WINDOW, dtype=np.int32)
    
    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom)

        sys.stderr.write("  reading midpoint counts\n")
        counts = mnase_track.get_nparray(chrom)

        sys.stderr.write("  smoothing counts with %dbp window\n" % 
                         MNASE_SMOOTH_WINDOW)
        
        weights = np.convolve(counts, window, mode='same')
                
        sys.stderr.write("  getting DNA sequence\n")
        dna = seq_track.get_seq_str(chrom)

        sys.stderr.write("  counting kmers\n")
        # count kmers on this chromosome
        genome.kmer.count_kmers(dna, kmer_size, weights, 
                                all_kmers, obs_kmers)


    seen = set([])
        
    sys.stdout.write("KMER REVCOMP.KMER MNASE.COUNT ALL.COUNT\n")
    for kmer_str in genome.kmer.list_kmers(kmer_size):
        if kmer_str in seen:
            continue
                
        if kmer_str in obs_kmers:
            obs_count = obs_kmers[kmer_str]
        else:
            obs_count = 0
            
        if kmer_str in all_kmers:
            all_count = all_kmers[kmer_str]
        else:
            obs_count = 0
            
        # combine counts with reverse complement kmer
        rc_kmer_str = genome.seq.revcomp(kmer_str)
        if rc_kmer_str != kmer_str:
            if rc_kmer_str in obs_kmers:
                obs_count += obs_kmers[rc_kmer_str]

            if rc_kmer_str in all_kmers:
                all_count += all_kmers[rc_kmer_str]
                
        sys.stdout.write("%s %s %d %d\n" % (kmer_str, rc_kmer_str, 
                                            obs_count, all_count))
        
        # flag that we have already seen this kmer
        # so that we don't report again when we see the revcomp kmer
        seen.add(kmer_str)
        seen.add(rc_kmer_str)

        

main()
            


        
    
