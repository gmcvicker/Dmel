import sys
import os

import genome.transcript
import genome.db
import genome.fastq

GENES_PATH = '/data/share/genes/dm3/flyBaseGene.txt'

HOME = os.environ['HOME']
FLYBASE_ID_FILE = HOME + "/data/Dmel/flybase/fbgn_annotation_ID_fb_2012_04.tsv.gz"

READ_LEN = 50

# name of one of the copies of H1, there are many identical copies...
H1_GENE_NAME = "CG31617-RA"


def get_h1_transcript(chrom_dict):
    sys.stderr.write("reading transcripts from %s\n" % GENES_PATH)
    trs = genome.transcript.read_transcripts(GENES_PATH, chrom_dict)

    for tr in trs:
        if tr.name == H1_GENE_NAME:
            return tr

    raise ValueError("could not find H1 gene\n")



def get_cdna(seq_track, tr):    
    cdna_str = ""

    for ex in tr.exons:
        ex_seq = seq_track.get_seq_str(ex.chrom, ex.start, ex.end)

        if tr.strand == -1:
            cdna_str += genome.seq.revcomp(ex_seq)
        else:
            cdna_str += ex_seq
        
    return cdna_str


def build_seq_dict(seq_str, read_len, allow_mismatch=True):
    start = 0
    end = start + read_len

    seq_dict = {}
    
    while end <= len(seq_str):
        s = seq_str[start:end]

        sys.stderr.write("%s (%d bp)\n" % (s, len(s)))

        if allow_mismatch:
            s_list = list(s)
            
            for i in range(read_len):
                orig_base = s_list[i]
                for new_base in ('A', "C", "G", "T"):
                    s_list[i] = new_base
                    s_mm = "".join(s_list)
                    seq_dict[s_mm] = True
                s_list[i] = orig_base

        seq_dict[s] = True
        start += 1
        end += 1

    return seq_dict


def main():

    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <fastq_file>\n")
        exit(2)

    fastq_file = sys.argv[1]
        
    gdb = genome.db.GenomeDB(assembly="dm3")

    seq_track = gdb.open_track("seq")

    chrom_dict = gdb.get_chromosome_dict()
    h1_tr = get_h1_transcript(chrom_dict)

    sys.stderr.write(">H1 (%s:%d-%d[%d])\n" % (h1_tr.chrom.name, h1_tr.start, h1_tr.end,
                                               h1_tr.strand))

    h1_seq = get_cdna(seq_track, h1_tr)

    sys.stderr.write("%s\n" % h1_seq)
    
    h1_dict = build_seq_dict(h1_seq, READ_LEN)
    
    match_count = 0
    dup_match_count = 0
    total_count = 0
    skipped_count = 0
    dup_match = {}
    
    for record in genome.fastq.read_fastq(fastq_file):
        read = record[1]
        if 'N' in read:
            skipped_count += 1
        else:
            total_count += 1

            if (total_count % 1000000) == 0:
                sys.stderr.write("\n%d\n" % total_count)

            if read in h1_dict:
                if read in dup_match:
                    dup_match_count += 1
                    sys.stderr.write("@")
                else:
                    dup_match[read] = True
                    match_count += 1
                    sys.stderr.write("#")
        
    sys.stdout.write("MATCH %d\n" % match_count)
    sys.stdout.write("DUP_MATCH: %d\n" % dup_match_count)
    sys.stdout.write("TOTAL %d\n" % total_count)
    
    seq_track.close()


main()
