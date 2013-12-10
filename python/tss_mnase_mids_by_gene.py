import sys
import argparse
import numpy as np


import math
import genome.db
import genome.transcript
import genome.gene

import gene_expr as ge


# TRACK_LIST = ['eye', 'haltere', 'leg', 'antenna', 'S2', 'S2_in_vitro']
TRACK_LIST = ['S2_in_vitro_031212', 'S2_in_vitro_080110']
# LINEAGES = ['dnase']

MNASE_TRACKS = {'eye' : 'mnase/mnase_midpoints_eye_122_to_159',
                'haltere' : 'mnase/mnase_midpoints_haltere_122_to_159',
                'leg' : 'mnase/mnase_midpoints_leg_122_to_159',
                'antenna' : 'mnase/mnase_midpoints_antenna_122_to_159',
                'S2' : 'mnase/mnase_midpoints_S2_101_to_191',
                'S2_in_vitro' : 'mnase/mnase_midpoints_S2_in_vitro_combined_101_to_191',
                'S2_in_vitro_031212' : 'mnase/mnase_midpoints_S2_in_vitro_031212_101_to_191',
                'S2_in_vitro_080110' : 'mnase/mnase_midpoints_S2_in_vitro_080110_101_to_191',
                'h1' : 'H1/h1_midpoints_101_to_191_combined',
                'hmgd' : 'HMGD/hmgd_midpoints_combined',
                'dnase' : "dnase/dnase_S2"}

FLANK = 1000


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("out_dir", help="directory to write output files to")

    args = parser.parse_args()
    
    return args

    

def float_str(val):
    if val is None:
        return "NA"    
    return "%.3f" % val
     

def write_gene_summary(gene_f, gene):
    gene_f.write("%s %s %d %d %d" %
                 (gene.name, gene.chrom.name,
                  gene.start, gene.end, gene.strand))

    for expr_sample in ge.EXPR_SAMPLES:
        gene_f.write(" %s %s" % (float_str(gene.expr[expr_sample]), 
                                 float_str(gene.expr_pct[expr_sample])))
    gene_f.write("\n")



def main():
    args = parse_args()
    
    gdb = genome.db.GenomeDB(assembly="dm3")
    
    chrom_dict = gdb.get_chromosome_dict()
    gene_dict = ge.get_genes(chrom_dict)

    n_sites = FLANK + FLANK + 1

    mnase_profiles = {}
    mnase_tracks = {}
    mnase_out_f = {}

    gene_summary_f = open(args.out_dir + "/" + "gene_summary.txt", "w")
    gene_summary_f.write("GENE.NAME CHROM START END STRAND")
    for expr_sample in ge.EXPR_SAMPLES:
        gene_summary_f.write(" EXPR.%s EXPR.PCT.%s" % (expr_sample, expr_sample))
    gene_summary_f.write("\n")
    
    for t in TRACK_LIST:
        mnase_profiles[t] = np.zeros(n_sites, np.uint32)
        mnase_tracks[t] = gdb.open_track(MNASE_TRACKS[t])
        out_path = args.out_dir + "/" + t + ".txt"
        mnase_out_f[t] = open(out_path, "w")
                    
    n_gene = 0
    
    for chrom in gdb.get_chromosomes():
        
        sys.stderr.write("%s\n" % chrom.name)
        genes = gene_dict[chrom.name]

        count = 0
        
        for gene in genes:            
            count += 1
            if count >= 10:
                sys.stderr.write(".")
                count = 0

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

            write_gene_summary(gene_summary_f, gene)
            
            for t in TRACK_LIST:
                track = mnase_tracks[t]
                vals = track.get_nparray(chrom.name, start, end)

                if gene.strand == -1:
                    vals = vals[::-1]

                out_f = mnase_out_f[t]
                out_f.write(" ".join([str(x) for x in vals]) + "\n")

                # write values for this gene
                mnase_profiles[t] += vals

    
    # close output files
    gene_summary_f.close()
    for f in mnase_out_f.values():
        f.close()
    
    for track in mnase_tracks.values():
        track.close()



main()    



