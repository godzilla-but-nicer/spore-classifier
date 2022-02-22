import numpy as np
import os
import sys
import time
from glob import glob

mode = sys.argv[1]
ethresh = float(sys.argv[2])

if mode == 'sample':
    prot_dir = 'data/sample_proteins'
elif mode == 'ww':
    prot_dir = 'data/labelled/ww_proteins'
db_path = 'data/reference_sequences/sporulation_sprot.fasta'
sample_faas = glob(prot_dir + '/*.faa')

print(f"Starting BLASTP on {len(sample_faas)} proteomes...")
overall_start = time.time()
for faa in sample_faas:
    # get the unique identifier for the MAG
    id = faa.split('/')[-1].split('.')[0]
    out_file = prot_dir + '/' + str(id) + '.txt'

    # make a directory for the output
    #if os.path.isdir(prot_dir + '/' + id):
    #    os.mkdir(prot_dir + '/' + id)

    # run the blast 
    start_time = time.time()
    print(f'\tBLASTP on proteome: {faa}')
    # best hit overhang from: https://www-ncbi-nlm-nih-gov.proxyiub.uits.iu.edu/books/NBK569839/
    os.system(f'blastp -query {faa} -db {db_path} -out {out_file} -num_threads 4 '
              f'-evalue {ethresh} -outfmt 6 -best_hit_overhang 0.25')
    print(f'\tCompleted in {(time.time() - start_time)} s')

overall_end = time.time()
print(f"Entire sample completed in {(overall_end - overall_start)/60} min")