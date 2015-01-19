

import os
import sys



DATA_DIR="/home/gmcvicker/data/Dmel/mod_encode/histone_modification/ChIP-chip/normalized-arrayfile_wiggle"

def main():
    filenames = os.listdir(DATA_DIR)

    num = 0

    seen = set([])
    
    for filename in filenames:
        # We only want the "Mvalues" wiggle files. According to
        # http://intermine.modencode.org/release-31/report.do?id=70000008
        # M values are The log-intensity ratio values (M-values) are
        # calculated for all perfect match (PM) probes as log2(ChIP
        # intensity) - log2(input intensity). The M values are then
        # shifted so that the mean is equal to 0.        
        if not filename.endswith(".Mvalues.wig.gz"):
            continue
        
        words = filename.split(":")

        mark_type = words[0]
        tissue = words[1]

        tissue_words = tissue.split("#")

        # We are only interested in experiments performed on S2 cells
        if not tissue_words[0].startswith("Cell-Line=S2"):
            continue

        cell_line = tissue_words[0].split("=")[1]
        cell_line = cell_line.replace("-", "_")

        replicate = words[3].replace("-", "_")

        postfix = words[-1]
        experiment_id = postfix.split(".")[0]

        key = "%s_%s_%s" % (mark_type, cell_line, experiment_id)
        if key in seen:
            sys.stderr.write("skipping duplicate experiment %s\n" % key)
            continue
        seen.add(key)

        num += 1
        sys.stdout.write("%d %s %s %s %s %s/%s\n" % (num, mark_type, cell_line, 
                                                     replicate, experiment_id,
                                                     DATA_DIR, filename))

        

main()
