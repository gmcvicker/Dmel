
import numpy as np
import genome.nuc as nuc

import sys

class DataTracks(object):
    def __init__(self, gdb, mnase_track):        
        self.mnase = gdb.open_track(mnase_track, "r")        
        self.seq = gdb.open_track("seq", "r")


    def close(self):
        self.mnase.close()
        self.seq.close()




class Profile(object):

    def __init__(self, length, midpoint, output_file):
        self.length = length
        self.mid = midpoint
        self.a = np.zeros(self.length, dtype=np.float32)
        self.c = np.zeros(self.length, dtype=np.float32)
        self.g = np.zeros(self.length, dtype=np.float32)
        self.t = np.zeros(self.length, dtype=np.float32)
        self.dinuc = np.zeros((self.length, nuc.N_NUC, nuc.N_NUC),
                              dtype=np.float32)
        
        self.mnase_mids = np.zeros(self.length, dtype=np.float32)

        self.out_f = output_file


    def add_dinuc_counts(self, seq_vals, weight=1.0):
        vals = np.zeros(seq_vals.size, dtype=np.int8)
        vals[seq_vals == ord('A')] = nuc.NUC_ID_A
        vals[seq_vals == ord('C')] = nuc.NUC_ID_C
        vals[seq_vals == ord('G')] = nuc.NUC_ID_G
        vals[seq_vals == ord('T')] = nuc.NUC_ID_T
        vals[seq_vals == ord('N')] = nuc.NUC_ID_N

        # combine first nucleotide, second nucleotide and position within
        # profile to create a flat index into the dicnucleotide count array
        idx1 = vals[:-1]
        idx2 = vals[1:]
        undef = (idx1 < 0) | (idx2 < 0)
        idx = (np.arange(vals.size-1) * nuc.N_NUC * nuc.N_NUC) + \
              (idx1 * nuc.N_NUC) + idx2
        self.dinuc.flat[idx[~undef]] += weight


    def add_counts(self, tracks, chrom, start, end, strand, weight=1.0,
                   max_mnase_val=20):
        mnase_vals = tracks.mnase.get_nparray(chrom, start, end)
        seq_vals = tracks.seq.get_nparray(chrom, start, end)

        # threshold positions where number of midpoints
        # exceeds pre-specified maximum. These may be amplification
        # or copy number artefacts
        if max_mnase_val:
            mnase_vals[mnase_vals > max_mnase_val] = max_mnase_val

        if strand == -1:
            # reverse complement data
            seq_vals = genome.seq.revcomp_nparray(seq_vals)
            mnase_vals = mnase_vals[::-1]

        is_a = seq_vals == ord('A')
        is_c = seq_vals == ord('C')
        is_g = seq_vals == ord('G')
        is_t = seq_vals == ord('T')

        self.a += is_a * weight
        self.c += is_c * weight
        self.g += is_g * weight
        self.t += is_t * weight

        self.add_dinuc_counts(seq_vals, weight=weight)

        self.mnase_mids += mnase_vals * weight




    def write(self, include_zero_pos=True):
        self.out_f.write("POS\tA\tC\tG\tT\tMNASE.MIDS")

        # write out dinucleotide portion of header
        for i in nuc.NUCS:
            for j in nuc.NUCS:
                self.out_f.write("\t%s%s" % (i, j))
        self.out_f.write("\n")
                
        for i in range(self.length):
            if include_zero_pos:
                pos = i - self.mid + 1
            else:
                # no 0th coordinate, skip from -1 to +1
                if i < self.mid:
                    pos = i - self.mid + 1
                else:
                    pos = i - self.mid + 2
            
            self.out_f.write("%d\t%g\t%g\t%g\t%g\t%g" %
                             (pos, self.a[i], self.c[i], self.g[i], self.t[i],
                              self.mnase_mids[i]))

            # write dinucleotide portion of counts
            for j in range(nuc.N_NUC):
                for k in range(nuc.N_NUC):
                    self.out_f.write("\t%g" % self.dinuc[i, j, k])

            self.out_f.write("\n")

