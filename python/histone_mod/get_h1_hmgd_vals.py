import sys
import genome.db
import genome.trackstat as trackstat
import numpy as np


REGION_SIZE = 2000

N_HMGD_REP = 5
N_H1_REP = 4

def get_h1_hmgd_tracks(gdb):    
    tracks = {}

    # open track for each HMGD replicate
    for i in range(N_HMGD_REP):
        track_name = "HMGD_rep%d" % (i+1)
        tracks[track_name] = gdb.open_track("HMGD/hmgd_midpoints_rep%d" % (i+1))

    # open track for each H1 replicate
    for i in range(N_H1_REP):
        track_name = "H1_rep%d" % (i+1)
        tracks[track_name] = gdb.open_track("H1/h1_midpoints_101_to_191_rep%d" %
                                            (i+1))

    # open track for MNase-seq, which is used as background
    tracks['mnase'] = gdb.open_track("mnase/mnase_midpoints_S2_101_to_191")
        
    return tracks
        
        

def get_track_scales(gdb, tracks):
    track_scales = {}
    
    for track_name in tracks.keys():
        track_stat = trackstat.get_stats(gdb, tracks[track_name])
        track_scales[track_name] = 1e9 / float(track_stat.sum)
        
        sys.stderr.write("%d mapped reads for track %s\n" %
                         (track_stat.sum, track_name))

    return track_scales

    

def main():
    gdb = genome.db.GenomeDB(assembly="dm3")

    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <file_with_coords>\n" % sys.argv[0])
        exit(2)
    
    tracks = get_h1_hmgd_tracks(gdb)

    track_scales = get_track_scales(gdb, tracks)

    track_names = tracks.keys()
    track_names.sort()

    input_file = open(sys.argv[1])

    vals = np.zeros(len(track_names), dtype=np.float32)
    
    # write header
    sys.stdout.write("CHROM START END %s\n" % " ".join(track_names))

    header = input_file.readline()
    
    for line in input_file:
        words = line.split()
        chrom_name = words[0]
        start = int(words[1])
        end = int(words[2])
        
        
        for i in range(len(track_names)):
            track = tracks[track_names[i]]
            track_scale = track_scales[track_names[i]]
            x = track.get_nparray(chrom_name, start, end)
            vals[i] = np.mean(x) * track_scale

        # only report regions where at least some values are defined
        sys.stdout.write("%s %d %d " % (chrom_name, start, end))
        sys.stdout.write(" ".join(["%.3f" % x for x in vals]))
        sys.stdout.write("\n")

    for track in tracks.values():
        track.close()
    
main()                
                               
            

            
            
        
        
        
        
    
    
