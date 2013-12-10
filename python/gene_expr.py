
import sys
import genome.gene
import genome.transcript

import math

GENES_PATH = "/data/share/genes/dm3/flyBaseGene.txt"



HIGH_EXPR_PCT = 0.75
LOW_EXPR_PCT = 0.25

EXPR_SAMPLES = ["CME_L1", "ML_DmD11", "ML_DmD17_c3", "ML_DmD20_c5", "S2_R+"]
EXPR_FILE = "/mnt/lustre/home/gmcvicker/data/Dmel/mod_encode/rna_seq/transcript_rna_seq_counts.txt"




def read_expr(sample):
    f = open(EXPR_FILE)

    header = f.readline()

    # get indices of RNA-seq counts, and group them by the
    # names of the expr samples
    words = header.split()
    idx = 0
    sample_idx = {}
    for col in words:
        if col.startswith("RNASEQ.RPKM"):
            tokens = col.split(".", 3)
            sample_name = tokens[2]

            sample_idx[sample_name] = idx
        idx += 1

    want_idx = sample_idx[sample]
    sys.stderr.write("using RNA-seq counts from sample: %s\n" % sample)

    expr_dict = {}

    for l in f:
        words = l.rstrip().split()
        tr_id = words[1]
        expr_dict[tr_id] = float(words[want_idx])

    f.close()

    return expr_dict



def set_gene_expr(genes):
    # create gene attributes for expr, keyed on sample name
    for g in genes:
        g.expr = {}
        g.expr_rank = {}
        g.expr_pct = {}
        g.has_high_expr = {}
        g.has_low_expr = {}

    # set gene expr for each sample
    for expr_sample in EXPR_SAMPLES:
        # get transcript expression for this sample
        tr_expr_dict = read_expr(expr_sample)

        for g in genes:
            # find max expr of all TFs
            max_expr = None
            g.name = g.transcripts[0].name

            for tr in g.transcripts:
                if tr.name in tr_expr_dict:
                    if max_expr is None or tr_expr_dict[tr.name] > max_expr:
                        max_expr = tr_expr_dict[tr.name]
                        g.name = tr.name
            g.expr[expr_sample] = max_expr

        
    # now rank genes based on their expression in each sample
    # in order to assign percentiles
    for expr_sample in EXPR_SAMPLES:
        ranked_genes = sorted(genes, key = lambda g: g.expr[expr_sample])

        n_ranked = 0
        for g in ranked_genes:
            if g.expr[expr_sample] is None:
                g.expr_rank[expr_sample] = 0
            else:
                n_ranked += 1
                g.expr_rank[expr_sample] = n_ranked

        # convert ranks to percentiles
        for g in ranked_genes:
            if g.expr_rank[expr_sample] == 0:
                g.expr_pct[expr_sample] = None
                g.has_high_expr[expr_sample] = None
                g.has_low_expr[expr_sample] = None
            else:
                g.expr_pct[expr_sample] = \
                  float(g.expr_rank[expr_sample]) / float(n_ranked)

                if g.expr_pct[expr_sample] > HIGH_EXPR_PCT:
                    g.has_high_expr[expr_sample] = True
                else:
                    g.has_high_expr[expr_sample] = False

                if g.expr_pct[expr_sample] <= LOW_EXPR_PCT:
                    g.has_low_expr[expr_sample] = True
                else:
                    g.has_low_expr[expr_sample] = False


def set_tr_expr(tr_list):
    """Sets expression attribute of each transcript in provided list
    of provided list of transcripts"""

    for tr in tr_list:
        tr.expr = {}
    
    # set gene expr for each sample
    for expr_sample in EXPR_SAMPLES:
        # get transcript expression for this sample
        expr_dict = read_expr(expr_sample)

        for tr in tr_list:
            if tr.name in expr_dict:
                tr.expr[expr_sample] = expr_dict[tr.name]
            else:
                tr.expr[expr_sample] = 0.0
                    

def get_transcripts(chrom_dict):
    """Returns dictionary of transcript lists, keyed on chromosome name.
    The returned transcripts have their expression attribute set"""
    
    sys.stderr.write("reading transcripts from %s\n" % GENES_PATH)
    tr_list = genome.transcript.read_transcripts(GENES_PATH, chrom_dict)
    set_tr_expr(tr_list)
    tr_dict = dict([(chrom_name, []) for chrom_name in chrom_dict.keys()])
    for tr in tr_list:
        tr_dict[tr.chrom.name].append(tr)

    return tr_dict

            
                    
def get_genes(chrom_dict):
    """Returns dictionary of transcript lists, keyed on chromosome name.
    The returned transcripts have their expression attribute set"""
    gene_dict = {}

    sys.stderr.write("reading transcripts from %s\n" % GENES_PATH)
    trs = genome.transcript.read_transcripts(GENES_PATH, chrom_dict)

    # sort transcripts and group them into genes
    sys.stderr.write("grouping transcripts into genes\n")
    genes = genome.gene.group_transcripts(trs)

    set_gene_expr(genes)

    for gene in genes:
        if gene.chrom.name in gene_dict:
            gene_dict[gene.chrom.name].append(gene)
        else:
            gene_dict[gene.chrom.name] = [gene]
    
    return gene_dict




