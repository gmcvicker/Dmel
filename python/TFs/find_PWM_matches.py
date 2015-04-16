
import sys
import numpy as np

import tables

import genome.db


PWM_FILENAME="/home/gmcvicker/data/Dmel/TFs/UmassPGFE_PWMfreq_PublicDatasetB_20150409.txt"

PSEUDOCOUNT = 1.0

N_NUC = 4
NUC_A = 0
NUC_C = 1
NUC_G = 2
NUC_T = 3


def read_n_lines(f, n):
    lines = []

    l = f.readline()
    
    while l:
        lines.append(l)

        if len(lines) == n:
            yield lines
            lines = []

        l = f.readline()

    return



def write_pwm(f, name, matrix):
    row_labels = ['A', 'C', 'G', 'T']

    n_row = len(row_labels)
    
    f.write(">%s\n" % name)
    if matrix.shape[0] != n_row:
        raise ValueError("expected %d rows\n" % n_row)

    for i in range(n_row):
        val_str = " ".join(["%.2f" % v for v in matrix[i,:]])
        f.write("%s | %s\n" % (row_labels[i], val_str))
    
    

def read_pwms():
    pwm_dict = {}

    # remember best info-content of PWM
    best_pwm_ic = {}

    f = open(PWM_FILENAME)
    
    # read blocks of 5 lines at a time
    block = read_n_lines(f, 5)

    for block in read_n_lines(f, 5):
        header = block[0]
        a_line = block[1]
        c_line = block[2]
        g_line = block[3]
        t_line = block[4]

        # sanity check for expected labels
        if not header.startswith(">"):
            raise ValueError("expected header to start with '>'")
        if not a_line.startswith("A |"):
            raise ValueError("expected line 1 to start with 'A |'")
        if not c_line.startswith("C |"):
            raise ValueError("expected line 2 to start with 'C |'")
        if not g_line.startswith("G |"):
            raise ValueError("expected line 3 to start with 'G |'")
        if not t_line.startswith("T |"):
            raise ValueError("expected line 4 to start with 'T |'")

        # get gene names from header
        header_tok = header[1:].split("_")
        if len(header_tok) < 3:
            raise ValueError("expected header to have 3 tokens")
        gene_name = header_tok[0].lower()
        fb_gene_name = header_tok[2]

        # parse out counts
        a_counts = [int(x) for x in a_line.split()[2:]]
        c_counts = [int(x) for x in c_line.split()[2:]]
        g_counts = [int(x) for x in g_line.split()[2:]]
        t_counts = [int(x) for x in t_line.split()[2:]]

        # make count matrix and add pseudocount
        count_matrix = np.array([a_counts, c_counts,
                                 g_counts, t_counts], dtype=np.float32)

        count_matrix += PSEUDOCOUNT

        # convert counts to frequencies
        count_totals = np.apply_along_axis(np.sum, 0, count_matrix)

        freq_matrix = count_matrix / count_totals
        
        # calc log ratio with background freqs
        # assume equal nucleotide freqs for background (0.25)
        pwm = np.log2(freq_matrix) - np.log2(0.25)

        # Information content, maximized when freqs are non-uniform
        info_content = np.sum(-freq_matrix * np.log2(freq_matrix))

        if gene_name in pwm_dict:
            # we already have a PWM for this gene...
            if info_content > best_pwm_ic[gene_name]:
                # information content of this PWM is better
                sys.stderr.write("replacing PWM for %s with one "
                                 "that has higher information content "
                                 "(%.2f > %.2f)\n" % 
                                 (gene_name, info_content,
                                  best_pwm_ic[gene_name]))

                pwm_dict[gene_name] = pwm
                best_pwm_ic[gene_name] = info_content
        else:
            # first time we've seen this gene
            pwm_dict[gene_name] = pwm
            best_pwm_ic[gene_name] = info_content

    return pwm_dict
    




def seq_to_base_matrix(seq_str):
    """creates a matrix with 1 row for each nucleotide, 1 column for each
    sequence position. Positions in sequence containing that nucleotide
    sequence containing that nucleotide are flagged with 1, other positions
    in that row are 0"""
    seq_vec = np.array(list(seq_str))
    seq_matrix = np.zeros((N_NUC, len(seq_str)), dtype=np.int8)

    seq_matrix[NUC_A, np.where(seq_vec == "A")[0]] = 1
    seq_matrix[NUC_C, np.where(seq_vec == "C")[0]] = 1
    seq_matrix[NUC_G, np.where(seq_vec == "G")[0]] = 1
    seq_matrix[NUC_T, np.where(seq_vec == "T")[0]] = 1

    return seq_matrix



def scan_pwm(seq_matrix, pwm):
    """calculate PWM scores across sequence represented by provided matrix.
    Returned scores are for when the PWM starts at that position and returned
    score matrix is of length len(sequence)-len(pwm)+1
    """

    seq_len = seq_matrix.shape[1]
    pwm_len = pwm.shape[1]

    scores = np.zeros(seq_len-pwm_len+1, dtype=np.float32)
    
    for i in range(N_NUC):
        # slide PWM over flags for each nucleotide
        # note that PWM needs to be reversed for convolution to work
        # correctly
        sys.stderr.write("    %d/%d\n" % (i+1, N_NUC))
        scores += np.convolve(seq_matrix[i,], pwm[i,::-1], mode='valid')
    
    return scores



def reverse_comp_pwm(pwm):
    """makes a reverse-complemented version of the PWM"""
    new_pwm = np.zeros(pwm.shape)

    new_pwm[NUC_A] = pwm[NUC_T, ::-1]
    new_pwm[NUC_C] = pwm[NUC_G, ::-1]
    new_pwm[NUC_G] = pwm[NUC_C, ::-1]
    new_pwm[NUC_T] = pwm[NUC_A, ::-1]

    return new_pwm



def create_carray(track, chrom):
    atom = tables.Float32Atom(dflt=0.0)        
    zlib_filter = tables.Filters(complevel=1, complib="zlib")
    
    # create CArray for this chromosome
    shape = [chrom.length]
    carray = track.h5f.createCArray(track.h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    return carray




def get_carray(h5f, chrom):
    return h5f.getNode("/%s" % chrom)








def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <TF_name>\n" %
                         sys.argv[0])

    gdb = genome.db.GenomeDB(assembly="dm3")
        
    tf_name = sys.argv[1].lower()
    
    pwm_dict = read_pwms()
    
    if tf_name not in pwm_dict:
        sys.stderr.write("No PWM for TF %s\n" % tf_name)
        exit(-1)
    
    pwm = pwm_dict[tf_name]
    rev_pwm = reverse_comp_pwm(pwm)
    sys.stderr.write("got PWM for TF %s:\n" % tf_name)
    write_pwm(sys.stderr, tf_name, pwm)
    sys.stderr.write("reverse complement\n")
    write_pwm(sys.stderr, tf_name, rev_pwm)
    
    fwd_track_name = "tf_pwm/" + tf_name + "_fwd"
    rev_track_name = "tf_pwm/" + tf_name + "_rev"
    fwd_pwm_track = gdb.create_track(fwd_track_name)
    rev_pwm_track = gdb.create_track(rev_track_name)
                    
    seq_track = gdb.open_track("seq")

    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom.name)

        fwd_carray = create_carray(fwd_pwm_track, chrom)
        rev_carray = create_carray(rev_pwm_track, chrom)

        sys.stderr.write("  getting sequence\n")
        seq_str = seq_track.get_seq_str(chrom)
        seq_matrix = seq_to_base_matrix(seq_str)
        
        sys.stderr.write("  computing PWM scores\n")
        pwm_scores = scan_pwm(seq_matrix, pwm)

        # scan reverse compliment PWM as well!
        sys.stderr.write("  computing reverse complement PWM scores\n")
        rev_pwm_scores = scan_pwm(seq_matrix, pwm)


        fwd_carray[0:pwm_scores.shape[0]] = pwm_scores
        rev_carray[0:rev_pwm_scores.shape[0]] = rev_pwm_scores
        
        # sys.stderr.write("  finding max score\n")
        # max_idx = np.argmax(pwm_scores)
        # sys.stderr.write("    max: %.2f, pos: %d\n" %
        #                  (pwm_scores[max_idx], max_idx+1))

        
    seq_track.close()
    fwd_pwm_track.close()
    rev_pwm_track.close()
    
    ####
    #### TODO: scan genome with PWM, save scores in HDF5
    #### then look at score distribution, decide on reasonable
    #### score threshold? or base score threshold on information 
    #### content?
    ####

    

main()
