import numpy as np
import pandas as pd
import sys
import os
from glob import glob


SAMPLE_SIZE = int(sys.argv[1])
MODE = sys.argv[2]
TAXONOMY = 'Firmicutes'

# set up partial paths depending on mode
if MODE == 'prot':
    big_path = 'data/gems_proteins'
    sample_path = 'data/sample_proteins'
    ext = 'faa'
elif MODE == 'nucl':
    big_path = 'data/gems_dna'
    sample_path = 'data/sample_dna'
    ext = 'nucl'

# load the metadata and use it to select taxa
meta_df = pd.read_csv('data/gems_metadata.tsv', sep='\t').dropna()
meta_taxon = meta_df[meta_df['ecosystem'].str.contains(TAXONOMY)]

# sample from the genome ids for that taxa
sample_genomes = np.random.choice(meta_taxon['genome_id'].values, size=SAMPLE_SIZE)

# manage directory
if not os.path.isdir(sample_path):
    os.mkdir(sample_path)
else:
    os.system(f'rm -r {sample_path}')
    os.mkdir(sample_path)

# perform copies
for genome in sample_genomes:
    os.system(f'cp {big_path}/{genome}.{ext} {sample_path}/')
