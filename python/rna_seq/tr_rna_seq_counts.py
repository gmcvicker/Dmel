import sys
import os

import gzip
import numpy as np

import genome.db
import genome.gene
import genome.transcript
import genome.trackstat




GENES_PATH = '/data/share/genes/dm3/flyBaseGene.txt'
# LINEAGES = ["CME_L1", "ML_DmD11", "ML_DmD17_c3", "ML_DmD20_c5",
#            "S2_R+"]

LABELS = ["S2_KD/Wt", "S2_KD/HMGD1KD", "S2_KD/H1KD"]


HOME = os.environ['HOME']
FLYBASE_ID_FILE = HOME + "/data/Dmel/flybase/fbgn_annotation_ID_fb_2012_04.tsv.gz"


def get_genes(chrom_dict):
    gene_dict = {}
    tr_dict = {}

    sys.stderr.write("reading transcripts from %s\n" % GENES_PATH)
    trs = genome.transcript.read_transcripts(GENES_PATH, chrom_dict)

    # sort transcripts and group them into genes
    sys.stderr.write("grouping transcripts into genes\n")
    genes = genome.gene.group_transcripts(trs)

    return genes



def read_gene_symbols():
    """reads associations between annotation IDs (CGs) and gene symbols"""

    f = gzip.open(FLYBASE_ID_FILE, "rb")

    id_dict = {}

    for line in f:
        line = line.rstrip()
        if line.startswith("#") or line == "":
            # skip blank and comment lines
            continue

        words = line.split("\t")
        symbol = words[0]
        anno_ids = words[3].split(",")

        if len(words) > 4:
            # there are 'secondary annotation IDs'
            anno_ids.extend(words[4].split(","))
            

        for anno_id in anno_ids:
            if anno_id in id_dict:
                id_dict[anno_id].append(symbol)
            else:
                id_dict[anno_id] = [symbol]
    f.close()

    return id_dict
    
    


def set_rna_seq_expr(tr, rna_seq_tracks, totals):
    """sets the read counts for each label on tr"""
    tr.read_counts = {}
    tr.expr = {}
    
    for label in LABELS:
        track = rna_seq_tracks[label]

        if not track.has_chromosome(tr.chrom):
            tr.read_counts[label] = None
            tr.expr[label] = None
            continue
        
        # count mapped reads that overlap exons
        total_read_count = 0

        for exon in tr.exons:
            vals = track.get_nparray(exon.chrom.name,
                                     exon.start, exon.end)
            total_read_count += np.sum(vals)

        tr.read_counts[label] = total_read_count

        # report expression as RPKM
        # reads per kilobase (of transcript) per million bases sequenced
        pseudocount = 1.0
        tr.expr[label] = ((total_read_count + pseudocount) * 1e9) /  \
            (float(tr.size() * totals[label]))
        


def choose_best_tr(gene):
    """Out of all the genes associated with this gene, chooses the one
    with the highest mean expression."""

    best_tr = None

    for tr in gene.transcripts:
        if (best_tr is None) or ((tr.expr is not None) and (tr.expr > best_tr.expr)):
            best_tr = tr

    return best_tr


def write_transcript(tr, totals, id_dict):
    gene_id = tr.name.split("-")[0]

    if gene_id in id_dict:
        gene_symbols = id_dict[gene_id]
        # there can be multiple symbols, just take first one
        gene_symbol = gene_symbols[0]
    else:
        # this gene has no symbol, just use gene id
        gene_symbol = gene_id
    
    # report transcript name and total length of exons
    sys.stdout.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t" % (gene_symbol, tr.name,
                                                       tr.chrom.name,
                                                       tr.start, tr.end,
                                                       tr.strand,
                                                       tr.size(),
                                                       tr.gc_count))
    
    for label in LABELS:
        if tr.read_counts[label] is None:
            sys.stdout.write("\tNA")
        else:
            sys.stdout.write("\t%g" % tr.read_counts[label])

    for label in LABELS:
        sys.stdout.write("\t%d" % totals[label])

    for label in LABELS:
        if tr.expr[label] is None:
            sys.stdout.write("\tNA")
        else:
            sys.stdout.write("\t%g" % tr.expr[label])
    

    sys.stdout.write("\n")



def get_total_counts(gdb, tracks):
    totals = {}
    
    for label in LABELS:
        stat = genome.trackstat.get_stats(gdb, tracks[label])
        total = stat.sum
        sys.stderr.write("%s total reads %d\n" % (label, total))
        totals[label] = total

    return totals



def set_gc_count(seq_track, tr):
    gc_count = 0
    
    for ex in tr.exons:
        seq_str = seq_track.get_seq_str(ex.chrom, ex.start, ex.end)
      
        for base in seq_str:
            if base in ("G", "C", "g", "c"):
                gc_count += 1
    tr.gc_count = gc_count
                
        
        
    

def main():
    gdb = genome.db.GenomeDB(assembly='dm3')
    
    chrom_dict = gdb.get_chromosome_dict()

    genes = get_genes(chrom_dict)

    rna_seq_tracks = {}

    seq_track = gdb.open_track("seq")

    for label in LABELS:
        track_name = "rna_seq/%s" % label
        rna_seq_tracks[label] = gdb.open_track(track_name)

    totals = get_total_counts(gdb, rna_seq_tracks)

    # write header
    headers = ["GENE.SYMBOL", "TRANSCRIPT.NAME", "CHROM", "START",
               "END", "STRAND", "TOTAL.EXON.LEN", "EXON.GC.COUNT"]

    for label in LABELS:
        lab = label.replace("/", ".")
        headers.append("RNASEQ.COUNT.%s" % lab)

    for label in LABELS:
        lab = label.replace("/", ".")
        headers.append("RNASEQ.TOTAL.MAPPED.%s" % lab)
        
    for label in LABELS:
        lab = label.replace("/", ".")
        headers.append("RNASEQ.RPKM.%s" % lab)

    sys.stdout.write("\t".join(headers) + "\n")

    id_dict = read_gene_symbols()

    for gene in genes:
        for tr in gene.transcripts:
            # count number of mapped RNA-seq reads overlapping exons
            set_rna_seq_expr(tr, rna_seq_tracks, totals)

        best_tr = choose_best_tr(gene)
        set_gc_count(seq_track, best_tr)
        write_transcript(best_tr, totals, id_dict)

    for track in rna_seq_tracks.values():
        track.close()        

    

main()
