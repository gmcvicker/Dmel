import sys

import numpy as np

import argparse

import genome.db
import genome.coord
import genome.transcript

GENES_PATH="/data/share/genes/dm3/flyBaseGene.txt"


def get_transcripts(chrom_dict):
    gene_dict = {}

    sys.stderr.write("reading transcripts from %s\n" % GENES_PATH)
    trs = genome.transcript.read_transcripts(GENES_PATH, chrom_dict)

    tr_dict = dict([(chrom_name, []) for chrom_name in chrom_dict.keys()])

    for tr in trs:
        tr_dict[tr.chrom.name].append(tr)

    return tr_dict



def get_histone_mod_tracks(gdb):
    f = open("/home/gmcvicker/data/Dmel/mod_encode/histone_modification/histone_mod_tasklist.txt")
    
    tracks = {}
    
    for line in f:
        words = line.rstrip().split()
        mark_name = words[1]
        cell_line = words[2]
        exp_id = words[4]

        track_name = "%s_%s_%s" % (mark_name, cell_line, exp_id)

        tracks[track_name] = gdb.open_track("histone_mods/" + track_name)

    f.close()
    
    return tracks

        
    


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--region_type", default="genomewide",
                        choices=("genomewide", "tss"))

    parser.add_argument("--region_size", default=2000, type=int)
    
    options = parser.parse_args()

    return options




def get_regions(gdb, options):
    regions = []
    
    if options.region_type == "genomewide":
        for chrom in gdb.get_chromosomes():
            start = 1
            end = start + options.region_size - 1
        
            while end <= chrom.length:
                c = genome.coord.Coord(chrom, start, end)
                regions.append(c)
                
                start = end + 1
                end = start + options.region_size - 1

    elif options.region_type == "tss":
        chrom_dict = gdb.get_chromosome_dict()
        tr_dict = get_transcripts(chrom_dict)

        seen_tss = set([])
        
        for chrom in gdb.get_chromosomes():
            tr_list = tr_dict[chrom.name]

            for tr in tr_list:
                if tr.strand == 1:
                    tss_pos = tr.start
                else:
                    tss_pos = tr.end

                # only want to count unique TSSs
                tss_name = "%s:%d" % (chrom.name, tss_pos)
                if tss_name in seen_tss:
                    continue
                seen_tss.add(tss_name)

                start = tss_pos - options.region_size / 2
                end = tss_pos + options.region_size / 2

                start = max(start, 1)
                end = min(chrom.length, end)
                
                c = genome.coord.Coord(chrom, start, end, strand=tr.strand)
                regions.append(c)

    return regions
        
    


def main():
    options = parse_args()
    
    gdb = genome.db.GenomeDB(assembly="dm3")
    
    tracks = get_histone_mod_tracks(gdb)

    track_names = tracks.keys()
    track_names.sort()
    
    # write header
    sys.stdout.write("CHROM START END N.DEF %s\n" % " ".join(track_names))

    vals = np.zeros(len(track_names), dtype=np.float32)
    counts = np.zeros(len(track_names), dtype=np.float32)

    region_list = get_regions(gdb, options)
    
    for region in region_list:
        n_def = None
        for i in range(len(track_names)):
            track = tracks[track_names[i]]
            x = track.get_nparray(region.chrom, region.start, region.end)

            isdef = ~(np.isnan(x))                
            if n_def is None:
                n_def = np.sum(isdef)

            if np.any(isdef):
                vals[i] = np.mean(x[isdef])
                # record minimum number of defined bases for 
                # all tracks for this region
                n_def = min(n_def, np.sum(isdef))
            else:
                vals[i] = np.nan

        # only report regions where at least some values are defined
        if n_def > 0:
            sys.stdout.write("%s %d %d %d " % 
                             (region.chrom.name, region.start, 
                              region.end, n_def))
            sys.stdout.write(" ".join(["%.3f" % x for x in vals]))
            sys.stdout.write("\n")

    for track in tracks.values():
        track.close()
    
main()                


            
            
        
        
        
        
    
    
