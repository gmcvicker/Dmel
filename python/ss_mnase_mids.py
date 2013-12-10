import sys
import argparse
import numpy as np
import gzip

import genome.db
import genome.transcript
import genome.gene
import genome.coord
import genome.seq

import splicesite
import gene_expr as ge


GENES_PATH = "/data/share/genes/hg19/ens_gene.txt"

EXPR_SAMPLE = "S2_R+"


MNASE_TRACKS = { 'mnase' : "mnase/mnase_midpoints_S2_101_to_191",
                 'hmgd' : 'HMGD/hmgd_midpoints_combined',
                 'h1' : 'H1/h1_midpoints_101_to_191_combined'}

FLANK = 200


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("out_dir", 
                        help="directory to write output files to")

    args = parser.parse_args()
    
    return args


def create_tss_mask(chrom, tr_list):
    mask = np.zeros(chrom.length, dtype=np.bool)

    # flag all TSSs. We want to be sure that 
    # we are looking at effects related to splicing, not from promoters.
    for tr in tr_list:
        if tr.strand == 1:
            s = max(tr.start - FLANK, 1)
            e = min(tr.start + FLANK, chrom.length)
        else:
            s = max(tr.end - FLANK, 1)
            e = min(tr.end + FLANK, chrom.length)

        mask[s-1:e] = 1

    return mask



def float_str(val):
    if val is None:
        return "NA"    
    return "%.3f" % val
     


def write_ss_vals(out_f, ss_list, val_track):
    for ss in ss_list:
        start = ss.start - FLANK
        end = ss.end + FLANK

        if start < 1 or end > ss.chrom.length:
            vals = np.empty(end - start + 1)
            vals[:] = np.nan
        else:
            vals = val_track.get_nparray(ss.chrom, start, end)
        
        if ss.strand == -1:
            vals = vals[::-1]

        if np.issubdtype(vals.dtype, float):
            # convert NaNs to 0
            vals = np.nan_to_num(vals)
        
        out_f.write(" ".join([str(x) for x in vals]) + "\n")



def write_ss_info_header(out_f):
    # TODO: want to add columns for IS.CANONICAL, NEAR.TSS
    out_f.write("CHROM POS STRAND TRANSCRIPT.ID TRANSCRIPT.EXPR "
                "SPLICE.SITE.TYPE EXON.TYPE NEAR.TSS DNA.SEQ IS.CANONICAL\n")

    
def write_ss_info(out_f, ss_list):
    for ss in ss_list:
        type_str = "NA"
        if ss.is_5_prime:
            type_str = "5p"
        elif ss.is_3_prime:
            type_str = "3p"
        
        out_f.write("%s %d %d %s %g %s %s %d %s %d\n" % 
                    (ss.chrom.name, ss.pos, ss.strand, ss.transcript_id,
                     ss.expr[EXPR_SAMPLE], type_str, ss.exon_type(), 
                     ss.near_tss, ss.seq_str, 
                     ss.is_canonical))
    


def set_ss_info(ss_list, seq_track, mask_tss): 
    """sets flag indicating whether splice site is close to a TSS, and
    also sets sequence of splice site and flags whether canonical or
    non-canonical"""

    # number of bases to get on each site of splice junction:
    n_bp = 4
    
    for ss in ss_list:
        # get sequence around splice site
        if ss.strand == 1:
            if ss.is_5_prime:
                start = ss.pos - n_bp
                end = ss.pos + n_bp - 1
            if ss.is_3_prime:
                start = ss.pos - n_bp + 1
                end = ss.pos + n_bp
                
            seq_str = seq_track.get_seq_str(ss.chrom, start, end)
                
        elif ss.strand == -1:
            if ss.is_5_prime:
                start = ss.pos - n_bp + 1
                end = ss.pos + n_bp
            if ss.is_3_prime:
                start = ss.pos - n_bp
                end = ss.pos + n_bp - 1
            
            seq_str = seq_track.get_seq_str(ss.chrom, start, end)
            seq_str = genome.seq.revcomp(seq_str)
                
        else:
            raise ValueError("unknown strand")

        # set exonic portion of sequence to uppercase, intronic
        # portion lowercase, and flag as canonical / non-canonical
        if ss.is_5_prime:
            exon_seq = seq_str[0:n_bp].upper()
            intron_seq = seq_str[n_bp:].lower()
            ss.seq_str = exon_seq + intron_seq
            ss.is_canonical = intron_seq.startswith("gt")
        elif ss.is_3_prime:
            intron_seq = seq_str[0:n_bp].lower()
            exon_seq = seq_str[n_bp:].upper()
            ss.seq_str = intron_seq + exon_seq
            ss.is_canonical = intron_seq.endswith("ag")
            
        # check if near to TSS
        start = ss.pos - FLANK
        end = ss.pos + FLANK
        ss.near_tss = np.any(mask_tss[start-1:end])
        
        


def main():
    args = parse_args()

    # open output files
    mnase_f = gzip.open(args.out_dir + "/mnase.txt.gz", "wb")
    hmgd_f = gzip.open(args.out_dir + "/hmgd.txt.gz", "wb")
    h1_f = gzip.open(args.out_dir + "/h1.txt.gz", "wb")
    info_f = gzip.open(args.out_dir + "/info.txt.gz", "wb")
    
    gdb = genome.db.GenomeDB(assembly="dm3")
    
    chrom_dict = gdb.get_chromosome_dict()
    
    tr_dict = ge.get_transcripts(chrom_dict)
    
    mnase_track = gdb.open_track(MNASE_TRACKS["mnase"])
    hmgd_track = gdb.open_track(MNASE_TRACKS['hmgd'])
    h1_track = gdb.open_track(MNASE_TRACKS['h1'])
    seq_track = gdb.open_track("seq")

    write_ss_info_header(info_f)

    chromosomes = gdb.get_chromosomes()
    #chromosomes = [gdb.get_chromosome('chr22')]
    
    for chrom in chromosomes:
        sys.stderr.write("%s\n" % chrom.name)
        tr_list = tr_dict[chrom.name]
        
        sys.stderr.write("creating TSS mask\n")
        mask_tss = create_tss_mask(chrom, tr_list)

        sys.stderr.write("  getting splice sites\n")
        ss_list = splicesite.create_splice_sites(tr_list)

        sys.stderr.write("  writing splice site info\n")
        set_ss_info(ss_list, seq_track, mask_tss)
        write_ss_info(info_f, ss_list)

        sys.stderr.write("  writing MNase counts\n")
        # write MNase-seq midpoint counts around splice sites
        write_ss_vals(mnase_f, ss_list, mnase_track)

        sys.stderr.write("  writing HMGD counts\n")
        # write MNase-seq midpoint counts around splice sites
        write_ss_vals(hmgd_f, ss_list, hmgd_track)

        sys.stderr.write("  writing H1 counts\n")
        # write MNase-seq midpoint counts around splice sites
        write_ss_vals(h1_f, ss_list, h1_track)
        
    seq_track.close()
    mnase_track.close()
    h1_track.close()
    hmgd_track.close()
    
    info_f.close()
    hmgd_f.close()
    h1_f.close()
    mnase_f.close()

    
main()    



