import csv
import os
import sys

HOME = os.environ['HOME']

TFBS_PATH = HOME + "/data/Dmel/redfly/redfly_tfbs.csv"


def get_tf_counts():
    """do a first pass through the file, counting number of unique occurances
    of each TFBS"""
    reader = csv.DictReader(open(TFBS_PATH, 'rb'))
    
    already_seen = set([])
    counts = {}
    
    for row in reader:
        tf_id = row['name']
        tf_coord = row['coordinates']
        
        key = tf_id  + ":" + tf_coord
        if key in already_seen:
            continue
        
        words = tf_coord.split(":")
        chrom_name = "chr" + words[0]
        (start_str, end_str) = words[1].split("..")

        if start_str == '' or end_str == '':
            continue

        start = int(start_str)
        end = int(end_str)

        if start > end:
            # bad coordinate
            continue

        words = tf_id.split(":")

        tf_name = words[0]
        
        already_seen.add(key)
        if tf_name in counts:
            counts[tf_name] += 1
        else:
            counts[tf_name] = 1
    
    return counts



def write_tfs(tf_counts):
    reader = csv.DictReader(open(TFBS_PATH, 'rb'))

    already_seen = set([])

    for row in reader:
        tf_id = row['name']
        tf_coord = row['coordinates']

        words = tf_coord.split(":")
        chrom_name = "chr" + words[0]
        (start_str, end_str) = words[1].split("..")

        if start_str == '' or end_str == '':
            continue

        start = int(start_str)
        end = int(end_str)

        if start > end:
            # skip this, bad coordinate
            continue
                        
        key = tf_id  + ":" + tf_coord
        
        if key in already_seen:
            continue

        already_seen.add(key)

        words = tf_id.split(":")
        tf_name = words[0]
        
        sys.stdout.write("%s\t%d\t%d\t%s\t%d\n" % 
                         (chrom_name, start, end, tf_name, tf_counts[tf_name]))
                
    

def main():
    tf_counts = get_tf_counts()
    
    write_tfs(tf_counts)
    


main()
