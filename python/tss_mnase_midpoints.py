import sys
import argparse
import numpy as np



import math
import genome.db
import genome.transcript
import genome.gene


# GENES_PATH="/data/share/genes/dm3/flyBaseGene.txt"
GENES_PATH = "/iblm/netapp/data1/external/Dmel/flyBaseGene.txt"

# LINEAGES = ['eye', 'haltere', 'leg', 'antenna', 'h1', 'hmgd']

# MNASE_TRACKS = {'eye' : 'mnase/mnase_midpoints_eye_122_to_159',
#                 'haltere' : 'mnase/mnase_midpoints_haltere_122_to_159',
#                 'leg' : 'mnase/mnase_midpoints_leg_122_to_159',
#                 'antenna' : 'mnase/mnase_midpoints_antenna_122_to_159',
#                 'h1' : 'H1/h1_midpoints_101_to_191_combined',
#                 'hmgd' : 'HMGD/hmgd_midpoints_combined'}

LINEAGES = ['eye', 'haltere', 'leg', 'antenna', 'S2_in_vitro_031212', 'S2']
MNASE_TRACKS = {'eye' : 'mnase/mnase_midpoints_eye_101_to_191',
                'haltere' : 'mnase/mnase_midpoints_haltere_101_to_191',
                'leg' : 'mnase/mnase_midpoints_leg_101_to_191',
                'antenna' : 'mnase/mnase_midpoints_antenna_101_to_191',
                'S2_in_vitro_031212' : 'mnase/mnase_midpoints_S2_in_vitro_031212_101_to_191',
                'S2' : 'mnase/mnase_midpoints_S2_101_to_191'}





FLANK = 1000

HIGH_EXPR_PCT = 0.75
LOW_EXPR_PCT = 0.25

# EXPR_FILE = "/mnt/lustre/home/gmcvicker/data/Dmel/mod_encode/rna_seq/transcript_rna_seq_counts.txt"


EXPR_FILE = "/iblm/netapp/data1/external/Dmel/mod_encode/transcript_rna_seq_counts.txt"





def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--high_expr", action="store_true",
                        default=False,
                        help="use high-expression genes only")

    parser.add_argument("--mid_expr", action="store_true",
                        default=False,
                        help="use mid-expression genes only")
    
    parser.add_argument("--low_expr", action="store_true",
                        default=False,
                        help="use low-expression genes only")

    args = parser.parse_args()

    if args.high_expr and args.low_expr:
        raise ValueError("cannot request both --low_expr and --high_expr")
    
    return args




def read_expr():

    sys.stderr.write("reading expression from %s\n" % EXPR_FILE)
    f = open(EXPR_FILE)

    header = f.readline()

    expr_dict = {}

    for l in f:
        words = l.rstrip().split()

        tr_id = words[0]

        #ttl_exon_len = float(words[6])
        #ttl_rna_seq_counts = sum(float(x) for x in words[7:])
        #expr = math.log((ttl_rna_seq_counts+1.0) / ttl_exon_len)

        rpkm = [float(x) for x in words[11:]]
        expr = sum(rpkm) / len(rpkm)

        expr_dict[tr_id] = expr

        sys.stderr.write("%s => %g\n" % (tr_id, expr))

    f.close()

    return expr_dict



def set_gene_tss_expr(genes):
    tr_expr_dict = read_expr()
    
    for gene in genes:
        gene.tss_expr = None
        
        for tr in gene.transcripts:
            # strip off trailing -RA -RB, etc
            tr_name = tr.name.split("-")[0]
            
            if tr_name in tr_expr_dict:
                if (tr.strand == 1 and tr.start == gene.start) or \
                    (tr.strand == -1 and tr.end == gene.end):

                    expr = tr_expr_dict[tr_name]

                    sys.stderr.write("EXPR: %s\n" % expr)
                    
                    # transcript that starts at TSS
                    if (gene.tss_expr is None) or (expr > gene.tss_expr):
                        gene.tss_expr = expr
            else:
                sys.stderr.write("%s not found in expr dict\n" % tr_name)
                        

    # now rank genes based on their expression
    genes = sorted(genes, key = lambda g: g.tss_expr)

    n_ranked = 0
    for g in genes:
        if g.tss_expr:
            n_ranked += 1
            g.tss_expr_rank = n_ranked
        else:
            g.tss_expr_rank = 0

    # convert ranks to percentiles
    for g in genes:
        if g.tss_expr_rank == 0:
            g.tss_expr_pct = None
            g.has_high_tss_expr = None
            g.has_low_tss_expr = None
        else:
            g.tss_expr_pct = g.tss_expr_rank / float(n_ranked)

            if g.tss_expr_pct > HIGH_EXPR_PCT:
                g.has_high_tss_expr = True
            else:
                g.has_high_tss_expr = False
            
            if g.tss_expr_pct <= LOW_EXPR_PCT:
                g.has_low_tss_expr = True
            else:
                g.has_low_tss_expr = False

            # sys.stderr.write("%g %g %s %s\n" % 
            #                  (g.tss_expr, g.tss_expr_pct,
            #                   str(g.has_high_tss_expr),
            #                   str(g.has_low_tss_expr)))

    

def get_genes(chrom_dict):
    gene_dict = {}

    sys.stderr.write("reading transcripts from %s\n" % GENES_PATH)
    trs = genome.transcript.read_transcripts(GENES_PATH, chrom_dict)

    # sort transcripts and group them into genes
    sys.stderr.write("grouping transcripts into genes\n")
    genes = genome.gene.group_transcripts(trs)

    set_gene_tss_expr(genes)

    for gene in genes:
        if gene.chrom.name in gene_dict:
            gene_dict[gene.chrom.name].append(gene)
        else:
            gene_dict[gene.chrom.name] = [gene]
    
    return gene_dict

    




def main():
    args = parse_args()
    
    gdb = genome.db.GenomeDB(assembly="dm3")
    
    chrom_dict = gdb.get_chromosome_dict()
    gene_dict = get_genes(chrom_dict)

    n_sites = FLANK + FLANK + 1

    mnase_profiles = {}
    mnase_tracks = {}
    for lineage in LINEAGES:
        mnase_profiles[lineage] = np.zeros(n_sites, np.uint32)
        mnase_tracks[lineage] = gdb.open_track(MNASE_TRACKS[lineage])

    if args.high_expr:
        sys.stderr.write("using high expression genes only\n")

    if args.low_expr:
        sys.stderr.write("using low expression genes only\n")

    if args.mid_expr:
        sys.stderr.write("using mid expression genes only\n")
    
    n_gene = 0
        
    for chrom in gdb.get_chromosomes():
        sys.stderr.write("%s\n" % chrom.name)
        genes = gene_dict[chrom.name]

        sys.stderr.write("there are %s genes\n" % len(genes))
        
        for gene in genes:
            if (gene.has_high_tss_expr is None) or \
                (gene.has_low_tss_expr is None):
                continue

            if args.high_expr and not gene.has_high_tss_expr:
                # skip genes that do not have high expression
                continue

            if args.low_expr and not gene.has_low_tss_expr:
                # only use low expression genes
                continue

            if args.mid_expr and (gene.has_high_tss_expr or 
                                  gene.has_low_tss_expr):
                continue
            
            if gene.strand == 1:
                tss = gene.start
            else:
                tss = gene.end
            
            start = tss - FLANK
            end = tss + FLANK

            if start < 1:
                continue
            if end > chrom.length:
                continue

            # keep track of the number of genes examined
            n_gene += 1

            for lineage in LINEAGES:
                track = mnase_tracks[lineage]
                vals = track.get_nparray(chrom.name, start, end)

                if gene.strand == -1:
                    vals = vals[::-1]

                mnase_profiles[lineage] += vals

    
    sys.stdout.write("POS\tN.SITE\t" + "\t".join(LINEAGES) + "\n")
    for i in range(n_sites):
        pos = i - FLANK

        sys.stdout.write("%d\t%d" % (pos, n_gene))
        for lineage in LINEAGES:
            sys.stdout.write("\t%d" % mnase_profiles[lineage][i])
        sys.stdout.write("\n")
                             
    
    for track in mnase_tracks.values():
        track.close()


    


    
main()    



