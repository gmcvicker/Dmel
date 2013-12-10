
import sys
import os

import gzip

HOME = os.environ['HOME']

FLYBASE_GENE_NAME_FILE = HOME + "/data/Dmel/flybase/fbgn_annotation_ID_fb_2012_04.tsv.gz"


def read_anno_id_to_flybase_gene_name():
    f = gzip.open(FLYBASE_GENE_NAME_FILE)

    anno2fb_dict = {}
    
    for line in f:
        line = line.rstrip()
        
        if line.startswith("#") or line == "":
            continue

        words = line.split("\t")

        if len(words) < 3:
            sys.stderr.write("WARNING: skipping line: %s\n" % line)
            continue
        
        fb_name = words[1]
        anno_id = words[3]

        if len(words) > 4:
            additional_ids = words[4].split()
        else:
            additional_ids = []

        anno_id_list = additional_ids + [anno_id]

        if fb_name:
            for a_id in anno_id_list:
                if a_id in anno2fb_dict:
                    sys.stderr.write("WARNING: repeated annotation id: %s\n" %
                                     a_id)
                else:
                    anno2fb_dict[a_id] = fb_name

                # sys.stderr.write("'%s' => '%s'\n" % (anno_id, fb_name))


    return anno2fb_dict
        


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <file_with_anno_ids_in_first_col.txt>\n" %
                         sys.argv[0])
        exit(2)
    
    anno2fb_dict = read_anno_id_to_flybase_gene_name()


    input_file = sys.argv[1]
    f = open(input_file)
    header = f.readline()
    if header.startswith("CG"):
        # doesn't look like a header
        f.rewind()
    else:
        sys.stdout.write(header)
    
    for line in f:
        words = line.rstrip().split(None, 1)
        anno_id = words[0]

        anno_id = anno_id.split("-")[0]
        
        if anno_id in anno2fb_dict:
            fb_name = anno2fb_dict[anno_id]
            sys.stdout.write("%s %s\n" % (fb_name, words[1]))
        else:
            sys.stdout.write(line)
            sys.stderr.write("could not find fbname for '%s'\n" % anno_id)
            
    



main()
