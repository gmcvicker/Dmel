
import sys
import tss
import genome.db

import numpy as np

TSS_FILE="/home/gmcvicker/data/Dmel/MNase/tss_mnase_midpoints_by_gene/gene_summary.fb_gene_names.txt"

# consensus TATA sequence
TATA_DEF = "TATAWAWR"

TATA_MAX_MISMATCH = 1

# look for TATA box contained in this region:
TATA_REGION_START = -50
TATA_REGION_END = -10


AMBI_CODES = {'H' : ('A', 'T', 'C'),
              'W' : ('A', 'T'),
              'R' : ('A', 'G'),
              'Y' : ('C', 'T'),
              'M' : ('A', 'C'),
              'K' : ('G', 'T'),
              'S' : ('G', 'C'),
              'N' : ('A', 'C', 'T', 'G')
              }



def find_tata_matches(tss_coord, seq_track):
    # get the TATA sequence 
    region = tss_coord.copy()
    if tss_coord.strand == 1:
        region.start = max(region.start + TATA_REGION_START - 1, 1)
        region.end = max(region.end + TATA_REGION_END - 1, 1)            
    else:
        region.start = min(region.start - TATA_REGION_END + 1, 
                           region.chrom.length)
        region.end = min(region.end - TATA_REGION_START + 1, 
                         region.chrom.length)
        
    region_seq = seq_track.get_seq_str(region.chrom, region.start, region.end)

    if tss_coord.strand == -1:
        # search negative strand sequence
        region_seq = genome.seq.revcomp(region_seq)

    match_list = []

    tata_len = len(TATA_DEF)

    # loop over region, looking for matches to TATA box
    for i in range(len(region_seq) - tata_len + 1):
        n_mismatch = 0

        for j in range(tata_len):
            # resolve ambiguity code in TATA definition
            nuc_code = TATA_DEF[j]
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

            if n_mismatch > TATA_MAX_MISMATCH:
                break
            
        if n_mismatch <= TATA_MAX_MISMATCH:
            # how far is end of matching string from TSS?
            if tss_coord.strand == 1:
                tss_dist = region.start - tss_coord.start + \
                    tata_len + i - 1

                match_end = tss_coord.start + tss_dist
                match_start = match_end - tata_len + 1
            else:
                # remember the match is to the reverse complement of the region
                tss_dist = tss_coord.start - (region.end - i - tata_len + 1)

                match_start = tss_coord.start - tss_dist
                match_end = match_start + tata_len - 1

            tata_match = genome.coord.Coord(tss_coord.chrom,
                                            match_start, match_end,
                                            strand=tss_coord.strand)

            tata_match.match = region_seq[i:i+tata_len]
            tata_match.n_mismatch = n_mismatch
            tata_match.tss_dist = tss_dist
            match_list.append(tata_match)

    return match_list
            
        
    
        


def write_header():
    sys.stdout.write("\t".join(["GENE.NAME", "CHROM", "TSS.CHROM.POS", 
                                "STRAND", "TATA.TSS.DIST", "TATA.CHROM.START",
                                "TATA.CHROM.END", "TATA.MATCH",
                                "TATAWAWR.CONSENSUS.DIST"]) + "\n")


        

def write_tata(f, tss_coord, tata):

    if tata is None:
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
                          tata.tss_dist, tata.start, tata.end,
                          tata.match, tata.n_mismatch))
    

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

    sys.stderr.write("getting list of TSSs\n")
    tss_list = read_tss_list(chrom_dict)
    
    seq_track = gdb.open_track("seq")

    write_header()
    
    for t in tss_list:
        # find matches to TATA sequence
        match_list = find_tata_matches(t, seq_track)

        # choose one with the fewest mismatches from canonical TATA sequence
        best_match = None
        for m in match_list:
            if (best_match is None) or (m.n_mismatch < best_match.n_mismatch):
                best_match = m

        
        write_tata(sys.stdout, t, best_match)
    
    seq_track.close()
        


main()

    
