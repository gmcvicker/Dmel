import sys
import genome.coord
import genome.db
import numpy as np
import genome.dist
import gene_expr as ge

DATA_TRACKS = {'mnase' : 'mnase/mnase_midpoints_S2_101_to_191',
                'h1' : 'H1/h1_midpoints_101_to_191_combined',
                'hmgd' : 'HMGD/hmgd_midpoints_combined',
                'dnase' : "dnase/dnase_S2"}

DHS_FILE = "/home/gmcvicker/data/Dmel/DNase/S2.dhs.high.txt"


FLANK = 1000

UNDEF_DIST = -2147483648

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
    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <output_dir>\n" % sys.argv[0])
        exit(2)
    
    out_dir = sys.argv[1]
    
    gdb = genome.db.GenomeDB(assembly="dm3")

    chrom_dict = gdb.get_chromosome_dict()

    dnase_regions = read_dnase_regions(chrom_dict)
    transcripts = ge.get_transcripts(chrom_dict)

    
    tracks = {}
    out_files = {}
    for track_label, track_name in DATA_TRACKS.items():
        tracks[track_label] = gdb.open_track(track_name)
        out_files[track_label] = open(out_dir + "/%s.txt" % track_label, "w")

    info_out_f = open(out_dir + "/dhs_info.txt", "w")    
    info_out_f.write("CHROM POS TSS.DIST\n")
    
    for chrom in gdb.get_chromosomes():
        # calculate distance to nearest TSS
        tr_list = transcripts[chrom.name]
        tss_dists = get_tss_dists(chrom, tr_list)
        
        for dhs in dnase_regions[chrom.name]:
            start = dhs.start - FLANK
            end = dhs.end + FLANK
            mid = (start + end) / 2

            if start < 1 or end > chrom.length:
                continue

            # write out information about DNase HS
            info_out_f.write("%s %d %d\n" % (chrom.name, mid, tss_dists[mid-1]))
            # write out MNase, ChIP-seq, DNase read counts
            for track_label in tracks.keys():
                vals = tracks[track_label].get_nparray(chrom, start, end)
                out_f = out_files[track_label]
                out_f.write(" ".join(["%d" % x for x in vals]) + "\n")
                


    info_out_f.close()
    for track_label in tracks.keys():
        tracks[track_label].close()
        out_files[track_label].close()
            



main()
