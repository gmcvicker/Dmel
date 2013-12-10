
import sys
import os

import numpy as np

import genome.db
import genome.transcript

TSS_REGION_SIZE = 250

GENES_PATH = "/data/share/genes/dm3/flyBaseGene.txt"

HOME = os.environ['HOME']
OUTPUT_DIR = HOME + "/data/Dmel/reads_by_region_type"

TRACKS = {'HMGD' : 'HMGD/hmgd_midpoints_combined',
#          'MNase_antenna' : 'mnase/mnase_midpoints_antenna_122_to_159',
#          'MNase_eye' : 'mnase/mnase_midpoints_eye_122_to_159',
#          'MNase_haltere' : 'mnase/mnase_midpoints_haltere_122_to_159',
#          'MNase_leg' : 'mnase/mnase_midpoints_leg_122_to_159',
          'MNase_S2' : 'mnase/mnase_midpoints_S2_101_to_191',
          'H1' : 'H1/h1_midpoints_101_to_191_combined'}

INTERGENIC = 0
EXON = 1
TSS_DOWNSTREAM = 2
TSS_UPSTREAM = 3
INTRON = 4
N_FLAGS = 5


def flag_to_type(flag):
    if flag == INTERGENIC:
        return "intergenic"
    elif flag == EXON:
        return "exon"
    elif flag == TSS_UPSTREAM:
        return "tss_upstream"
    elif flag == TSS_DOWNSTREAM:
        return "tss_downstream"
    elif flag == INTRON:
        return "intron"
    


def get_transcripts(chrom_dict):
    tr_dict = dict([(x, []) for x in chrom_dict.keys()])
    
    tr_list = genome.transcript.read_transcripts(GENES_PATH, chrom_dict)

    # group transcripts by chromosome
    for tr in tr_list:
        tr_dict[tr.chrom.name].append(tr)

    return tr_dict


    
def flag_annotations(chrom, tr_list):

    flags = np.zeros(chrom.length)

    flags[:] = INTERGENIC

    # flag intronic sequences
    for tr in tr_list:
        introns = tr.get_introns()

        for intron in tr.get_introns():
            flags[intron.start-1:intron.end] = INTRON

    # flag exons
    for tr in tr_list:
        for exon in tr.exons:
            flags[exon.start-1:exon.end] =  EXON

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

        flags[tss_ds_start-1:tss_ds_end] = TSS_DOWNSTREAM
        flags[tss_us_start-1:tss_us_end] = TSS_UPSTREAM        

    return flags
    


def main():
    gdb = genome.db.GenomeDB(assembly="dm3")

    # one output file for each input track
    track_dict = {}
    output_files = {}
    total_sites = {}
    total_counts = {}

    for track_name, track_path in TRACKS.items():
        track_dict[track_name] = gdb.open_track(track_path)

        output_files[track_name] = \
          open("%s/%s.txt" % (OUTPUT_DIR, track_name), "w")

        total_sites[track_name] = dict([(x , 0) for x in range(N_FLAGS)])
        total_counts[track_name] = dict([(x , 0) for x in range(N_FLAGS)])

    
    sys.stderr.write("reading transcripts\n")
    chrom_dict = gdb.get_chromosome_dict()
    tr_dict = get_transcripts(chrom_dict)
    
    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom.name)
        sys.stderr.write("  flagging annotations\n")
        flags = flag_annotations(chrom, tr_dict[chrom.name])

        for track_name in track_dict.keys():
            track = track_dict[track_name]
            read_counts = track.get_nparray(chrom).astype(np.int32)

            sys.stderr.write("  counting reads of type %s\n" % track_name)
            for flag in range(N_FLAGS):
                idx = np.where(flags == flag)[0]
                
                total_sites[track_name][flag] += idx.size
                total_counts[track_name][flag] += np.sum(read_counts[idx])


    # write output
    for track_name in track_dict.keys():
        f = output_files[track_name]
        f.write("REGION.TYPE TOTAL.SITES READ.COUNT\n")

        for flag in range(N_FLAGS):
            type_name = flag_to_type(flag)
            f.write("%s %d %d\n" % (type_name, total_sites[track_name][flag],
                                    total_counts[track_name][flag]))

    for track in track_dict.values():
        track.close()

    for f in output_files.values():
        f.close()
                                         
    
            

main()
        
        
