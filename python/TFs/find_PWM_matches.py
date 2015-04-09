
import sys
import numpy as np

PWM_FILENAME="/home/gmcvicker/data/Dmel/TFs/UmassPGFE_PWMfreq_PublicDatasetB_20150409.txt"

PSEUDOCOUNT = 1.0


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
    



def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <TF_name>\n" %
                         sys.argv[0])
    
    tf_name = sys.argv[1].lower()
    
    pwm_dict = read_pwms()
    
    if tf_name not in pwm_dict:
        sys.stderr.write("No PWM for TF %s\n" % tf_name)
        exit(-1)
    
    pwm = pwm_dict[tf_name]
    sys.stderr.write("got PWM for TF %s:\n" % tf_name)
    write_pwm(sys.stderr, tf_name, pwm)

    ####
    #### TODO: scan genome with PWM, save scores in HDF5
    #### then look at score distribution, decide on reasonable
    #### score threshold? or base score threshold on information 
    #### content?
    ####


    

    

main()
