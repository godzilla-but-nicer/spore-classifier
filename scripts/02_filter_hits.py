import pandas as pd
import os
import sys
import time
from glob import glob

mode = sys.argv[1]

if mode == 'sample':
    prot_dir = 'data/sample_proteins'
elif mode =='ww':
    prot_dir = 'data/labelled/ww_proteins'
elif mode == 'browne':
    prot_dir = 'data/labelled/browne_proteins'
    

blast_outputs = glob(prot_dir + '/*.txt')
for bof_path in blast_outputs:
    # for each file we're going to collect a list of output lines
    filtered_output = []
    with open(bof_path, 'r') as bof:
        for li, line in enumerate(bof.readlines()):
            spline = line.split('\t')
            # initialize the variables for the search
            if li == 0:
                query_seq = spline[0]  # seq label in query file
                best_escore = float(spline[10])
                best_hit = line
            # if we find another match with better score
            elif query_seq == spline[0] and best_escore > float(spline[10]):
                best_hit = line
                best_escore = float(spline[10])
            # if were looking at a different sequence add to list and reinitialize
            elif query_seq != spline[0]:
                filtered_output.append(best_hit)
                query_seq = spline[0]  # seq label in query file
                best_escore = float(spline[10])
                best_hit = line
    
    # now we want to write a new file
    if not os.path.isdir(prot_dir + '/best_hits'):
        os.mkdir(prot_dir + '/best_hits')
    genome_id = bof_path.split('/')[-1].split('.txt')[0]
    outfile_path = f'{prot_dir}/best_hits/{genome_id}_best.txt'
    with open(outfile_path, 'w') as outf:
        outf.writelines(filtered_output)