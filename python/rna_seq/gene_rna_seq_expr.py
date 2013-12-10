
import genome.db

import numpy as np

import genome.gene
import genome.transcript
import sys
import os

HOME = os.environ['HOME']

GENES_PATH = '/data/share/genes/dm3/flyBaseGene.txt'

TR_RNA_SEQ_PATH = HOME + "/data/Dmel/mod_encode/rna_seq/transcript_rna_seq_counts.txt"

LINEAGES = ["CME_L1", "ML_DmD11", "ML_DmD17_c3", "ML_DmD20_c5", "S2_R+"]





def main():
    gdb = genome.db.GenomeDB(assembly='dm3')

    f = open(TR_RNA_SEQ_PATH, "r")

    header = f.readline()

    gene_best_tr = {}
    gene_highest_expr = {}

    for line in f:
        words = line.rstrip().split()

        tr_name = words[0]
        gene_name = tr_name.split("-")[0]

        tr_len = int(words[5])

        rna_seq_counts = np.array([int(x) for x in words[6:]])

        expr = float(np.sum(rna_seq_counts)) / tr_len

        sys.stderr.write("mean_expr: %g\n" % expr)

        if gene_name in gene_best_tr:
            max_expr = gene_highest_expr[gene_name]
        else:
            max_expr = None

        if max_expr is None or expr > max_expr:
            # this transcript has highest mean expression
            # for this gene
            gene_highest_expr[gene_name] = expr
            gene_best_tr[gene_name] = tr_name

    
    # do second pass through file only reporting transcript with highest
    # expression (avg'd across cell lines)
    f.seek(0)
    header = f.readline()

    for line in f:
        words = line.split()
        tr_name = words[0]

        gene_name = tr_name.split("-")[0]

        if tr_name == gene_best_tr[gene_name]:
            sys.stderr.write(line)

    f.close()
    
                
                
        
        

    

main()
