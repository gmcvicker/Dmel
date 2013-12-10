
import sys
import tss
import genome.db

import argparse

import numpy as np

TSS_FILE="/home/gmcvicker/data/Dmel/MNase/tss_mnase_midpoints_by_gene/gene_summary.fb_gene_names.txt"



AMBI_CODES = {'H' : ('A', 'T', 'C'),
              'V' : ('A', 'C', 'G'),
              'D' : ('A', 'G', 'T'),
              'B' : ('C', 'G', 'T'),
              'W' : ('A', 'T'),
              'R' : ('A', 'G'),
              'Y' : ('C', 'T'),
              'M' : ('A', 'C'),
              'K' : ('G', 'T'),
              'S' : ('G', 'C'),
              'N' : ('A', 'C', 'T', 'G')
              }



def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("region_start",  type=int,
                        help="start of region to scan (relative to TSS)")
    
    parser.add_argument("region_end", type=int, 
                        help="end of region to scan (relative to TSS)")
    
    parser.add_argument("--max_mismatch", type=int, default=1,
                        help="max number of motif mismatches")

    parser.add_argument("motif", help="motif to scan for (e.g. TATAWAWR)")

    args = parser.parse_args()

    if args.region_start > args.region_end:
        raise ValueError("region start must be <= region end")
    
    return args
    

    


def find_motif_matches(tss_coord, seq_track, region_start, region_end, 
                      motif, max_mismatch):
    # get the sequence around the TSS
    region = tss_coord.copy()
    if tss_coord.strand == 1:
        region.start = max(region.start + region_start - 1, 1)
        region.end = max(region.end + region_end - 1, 1)            
    else:
        region.start = min(region.start - region_end + 1, 
                           region.chrom.length)
        region.end = min(region.end - region_start + 1, 
                         region.chrom.length)
        
    region_seq = seq_track.get_seq_str(region.chrom, region.start, 
                                       region.end)

    if tss_coord.strand == -1:
        # search negative strand sequence
        region_seq = genome.seq.revcomp(region_seq)

    match_list = []

    motif_len = len(motif)

    # loop over region, looking for matches to motif
    for i in range(len(region_seq) - motif_len + 1):
        n_mismatch = 0

        for j in range(motif_len):
            # resolve ambiguity codes in motif definition
            nuc_code = motif[j]
            if nuc_code in AMBI_CODES:
                nucs = AMBI_CODES[nuc_code]
            else:
                nucs = (nuc_code)

            # count mismatches
            is_match = False
            for nuc in nucs:
                if nuc == region_seq[i+j]:
                    is_match = True
                    break
            if not is_match:
                n_mismatch += 1

            if n_mismatch > max_mismatch:
                break
            
        if n_mismatch <= max_mismatch:
            # how far is end of matching string from TSS?
            if tss_coord.strand == 1:
                tss_dist = region.start - tss_coord.start + \
                    motif_len + i - 1

                match_end = tss_coord.start + tss_dist
                match_start = match_end - motif_len + 1
            else:
                # remember the match is to the reverse 
                # complement of the region
                tss_dist = tss_coord.start - (region.end - i - motif_len + 1)

                match_start = tss_coord.start - tss_dist
                match_end = match_start + motif_len - 1

            match = genome.coord.Coord(tss_coord.chrom,
                                       match_start, match_end,
                                       strand=tss_coord.strand)

            match.seq = region_seq[i:i+motif_len]
            match.n_mismatch = n_mismatch
            match.tss_dist = tss_dist
            match_list.append(match)

    return match_list
            
        
    
        


def write_header(args):
    sys.stdout.write("# region_start: %d\n" % args.region_start)
    sys.stdout.write("# region_end: %d\n" % args.region_end)
    sys.stdout.write("# max_mismatch: %d\n" % args.max_mismatch)
    sys.stdout.write("# motif: %s\n" % args.motif)
    
    sys.stdout.write("\t".join(["GENE.NAME", "CHROM", "TSS.CHROM.POS", 
                                "STRAND", "MOTIF.TSS.DIST", 
                                "MOTIF.CHROM.START",
                                "MOTIF.CHROM.END", "MOTIF.SEQ",
                                "MOTIF.MISMATCH"]) + "\n")


        

def write_match(f, tss_coord, match):

    if match is None:
        sys.stdout.write("%s\t%s\t%d\t%d\t"
                         "NA\tNA\tNA\t"
                         "NA\tNA\n" % \
                         (tss_coord.name, tss_coord.chrom.name, 
                          tss_coord.start, tss_coord.strand))
    else:
        sys.stdout.write("%s\t%s\t%d\t%d\t"
                         "%d\t%d\t%d\t"
                         "%s\t%d\n" % \
                         (tss_coord.name, 
                          tss_coord.chrom.name, tss_coord.start, 
                          tss_coord.strand,
                          match.tss_dist, match.start, match.end,
                          match.seq, match.n_mismatch))
    

def read_tss_list(chrom_dict):
    f = open(TSS_FILE, "r")
    header = f.readline()

    tss_list = []
    
    for line in f:
        words = line.split()

        name = words[0]
        chrom = chrom_dict[words[1]]
        gene_start = int(words[2])
        gene_end = int(words[3])
        strand = int(words[4])

        if strand == 1:
            tss_pos = gene_start
        elif strand == -1:
            tss_pos = gene_end
        else:
            raise ValueError("unknown strand")
        
        t = genome.coord.Coord(chrom, tss_pos, tss_pos, strand=strand)
        t.name = name
        
        tss_list.append(t)
        
    f.close()

    return tss_list
    
    


def main():    
    gdb = genome.db.GenomeDB(assembly="dm3")

    chrom_dict = gdb.get_chromosome_dict()

    args = parse_args()
    
    sys.stderr.write("getting list of TSSs\n")
    tss_list = read_tss_list(chrom_dict)
    
    seq_track = gdb.open_track("seq")

    write_header(args)
    
    for t in tss_list:
        # find matches to motif sequence
        match_list = find_motif_matches(t, seq_track,
                                        args.region_start,
                                        args.region_end,
                                        args.motif,
                                        args.max_mismatch)
                                        

        # choose one with the fewest mismatches to motif
        best_match = None
        for m in match_list:
            if (best_match is None) or (m.n_mismatch < best_match.n_mismatch):
                best_match = m

        
        write_match(sys.stdout, t, best_match)
    
    seq_track.close()
        


main()

    
