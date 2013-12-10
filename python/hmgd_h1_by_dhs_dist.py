import sys
import genome.coord
import genome.db
import gene_expr as ge
import numpy as np
import genome.dist

DATA_TRACKS = {'mnase' : 'mnase/mnase_midpoints_S2_101_to_191',
                'h1' : 'H1/h1_midpoints_101_to_191_combined',
                'hmgd' : 'HMGD/hmgd_midpoints_combined',
                'dnase' : "dnase/dnase_S2"}

DHS_FILE = "/home/gmcvicker/data/Dmel/DNase/S2.dhs.high.txt"


TSS_REGION_SIZE = 250
WIN_SZ = 1000

UNDEF_DIST = -2147483648


FLAG_INTERGENIC = 0
FLAG_EXON = 1
FLAG_TSS_DOWNSTREAM = 2
FLAG_TSS_UPSTREAM = 3
FLAG_INTRON = 4
N_FLAGS = 5

FLAG_NAMES = ['intergenic', 'exon', 'tss_upstream', 'tss_downstream',
              'intron']


def flag_annotations(chrom, tr_list):

    flags = np.zeros(chrom.length)

    flags[:] = FLAG_INTERGENIC

    # flag intronic sequences
    for tr in tr_list:
        introns = tr.get_introns()

        for intron in tr.get_introns():
            flags[intron.start-1:intron.end] = FLAG_INTRON

    # flag exons
    for tr in tr_list:
        for exon in tr.exons:
            flags[exon.start-1:exon.end] =  FLAG_EXON

    # flag promoter region
    for tr in tr_list:
        if tr.strand == 1:
            # divide promoter region into upstream and downstream
            # of TSS
            tss_us_start = tr.start - TSS_REGION_SIZE
            tss_us_end = tr.start - 1

            tss_ds_start = tr.start
            tss_ds_end = tr.start + TSS_REGION_SIZE - 1
        elif tr.strand == -1:
            tss_us_start = tr.start + 1
            tss_us_end = tr.start + TSS_REGION_SIZE

            tss_ds_start = tr.start - TSS_REGION_SIZE + 1
            tss_ds_end = tr.start            
        else:
            raise ValueError("unknown strand")

        flags[tss_ds_start-1:tss_ds_end] = FLAG_TSS_DOWNSTREAM
        flags[tss_us_start-1:tss_us_end] = FLAG_TSS_UPSTREAM        

    return flags

    



def read_dnase_regions(chrom_dict):
    regions = dict([(x, []) for x in chrom_dict.keys()])

    f = open(DHS_FILE, "r")
    header = f.readline()
    
    for l in f:
        words = l.split()
        chrom = chrom_dict["chr" + words[0]]
        pos = int(words[1])

        c = genome.coord.Coord(chrom, pos, pos)
        regions[chrom.name].append(c)

    return regions

    
    

def get_dhs_dists(chrom, dhs_list):
    dhs_dists = np.empty(chrom.length, dtype=np.int32)
    dhs_dists[:] = UNDEF_DIST

    # flag location of transcript start
    for dhs in dhs_list:
        dhs_dists[dhs.start-1:dhs.end] = 0

    return genome.dist.calc_dists_array(dhs_dists)



def get_tss_dists(chrom, tr_list):
    tss_dists = np.empty(chrom.length, dtype=np.int32)
    tss_dists[:] = UNDEF_DIST

    # flag location of transcript start
    for tr in tr_list:
        if tr.strand == 1:
            tss_dists[tr.start-1] = 0
        else:
            tss_dists[tr.end-1] = 0

    return genome.dist.calc_dists_array(tss_dists)



def main(): 
    gdb = genome.db.GenomeDB(assembly="dm3")

    chrom_dict = gdb.get_chromosome_dict()

    dnase_regions = read_dnase_regions(chrom_dict)
    transcripts = ge.get_transcripts(chrom_dict)
    
    tracks = {}
    for track_label, track_name in DATA_TRACKS.items():
        tracks[track_label] = gdb.open_track(track_name)

    # write header
    sys.stdout.write("CHROM START END MEAN.DHS.DIST MIN.DHS.DIST "
                     "MEAN.TSS.DIST MIN.TSS.DIST ")
    sys.stdout.write(" ".join(["%s.COUNT" % x.upper() for x in tracks.keys()]))
    sys.stdout.write(" ")
    sys.stdout.write(" ".join(["%s.COUNT" % x.upper() for x in FLAG_NAMES]))
    sys.stdout.write("\n")
    
    
    for chrom in gdb.get_chromosomes():
        # calculate distance to nearest DHS
        sys.stderr.write("%s\n" % chrom.name)

        if len(dnase_regions) == 0:
            # no DHSs on this chromosome
            continue

        sys.stderr.write("  flagging annotations\n")
        flags = flag_annotations(chrom, transcripts[chrom.name])
        
        dhs_dists = get_dhs_dists(chrom, dnase_regions[chrom.name])
        tss_dists = get_tss_dists(chrom, transcripts[chrom.name])

        for win_start in range(1, chrom.length, WIN_SZ):
            win_end = win_start + WIN_SZ-1

            if win_end > chrom.length:
                continue

            mean_dhs_dist = np.mean(dhs_dists[win_start-1:win_end].astype(np.float32))
            min_dhs_dist = np.min(dhs_dists[win_start-1:win_end])

            
            mean_tss_dist = np.mean(tss_dists[win_start-1:win_end].astype(np.float32))
            min_tss_dist = np.min(tss_dists[win_start-1:win_end])

            # write out information about region
            sys.stdout.write("%s %d %d %.2f %d %.2f %d" % 
                             (chrom.name, win_start, win_end, mean_dhs_dist, 
                              min_dhs_dist, mean_tss_dist, min_tss_dist))
            
            # write out MNase, ChIP-seq, DNase read counts
            for track_label in tracks.keys():
                vals = tracks[track_label].get_nparray(chrom, win_start, win_end)
                sys.stdout.write(" %d" % np.sum(vals))

            # write out number of bases of each annotation type
            for flag_id in range(N_FLAGS):
                n_flagged = np.sum((flags[win_start-1:win_end] == flag_id))
                sys.stdout.write(" %d" % n_flagged)
                
            sys.stdout.write("\n")
                
    for track_label in tracks.keys():
        tracks[track_label].close()

            



main()
