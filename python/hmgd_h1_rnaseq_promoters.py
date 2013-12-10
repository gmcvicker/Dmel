import sys
import genome.coord
import genome.db
import gene_expr as ge
import numpy as np
import genome.dist

DATA_TRACKS = {'mnase' : 'mnase/mnase_midpoints_S2_101_to_191',
                'h1'   : 'H1/h1_midpoints_101_to_191_combined',
                'hmgd' : 'HMGD/hmgd_midpoints_combined'}


INPUT_FILE = "/home/gmcvicker/data/Dmel/KD_RNA_seq/tr_rna_seq_counts.txt"
    
TSS_FLANKING = 500

def read_tss_regions(chrom_dict):
    regions = []

    f = open(INPUT_FILE, "r")
    header = f.readline()
    
    for l in f:
        words = l.split()

        # read transcript coordinates

        chrom = chrom_dict[words[2]]
        tr_strand = genome.coord.parse_strand(words[5])
        tr_start = int(words[3])
        tr_end = int(words[4])
        
        if tr_strand == 1:
            start = tr_start - TSS_FLANKING
            end = tr_start + TSS_FLANKING
        else:
            start = tr_end - TSS_FLANKING
            end = tr_end + TSS_FLANKING

        if start < 1:
            start = 1
        if end > chrom.length:
            end = chrom.length
            
        region = genome.coord.Coord(chrom, start, end)
        regions.append(region)

        region.gene_symbol = words[0]
        

    return regions

    
    



def main(): 
    gdb = genome.db.GenomeDB(assembly="dm3")

    chrom_dict = gdb.get_chromosome_dict()

    tss_regions = read_tss_regions(chrom_dict)
    
    tracks = {}
    for track_label, track_name in DATA_TRACKS.items():
        tracks[track_label] = gdb.open_track(track_name)

    # write header
    sys.stdout.write("CHROM START END GENE.SYMBOL ")
    sys.stdout.write(" ".join(["%s.COUNT" % x.upper() for x in tracks.keys()]))
    sys.stdout.write("\n")
        

    for region in tss_regions:
        # write out information about region
        sys.stdout.write("%s %d %d %s" % 
                         (region.chrom.name, region.start, region.end,
                          region.gene_symbol))

        # write out MNase, ChIP-seq read counts
        for track_label in tracks.keys():
            track = tracks[track_label]
            if track.has_chromosome(region.chrom):
                vals = track.get_nparray(region.chrom, region.start, region.end)
                sys.stdout.write(" %d" % np.sum(vals))
            else:
                sys.stdout.write(" NA")

        sys.stdout.write("\n")
                
    for track_label in tracks.keys():
        tracks[track_label].close()

            



main()
