#!/usr/bin/env python

import pandas as pd
import re
import sys
from glob import glob
from Bio import SeqIO

hit_paths = glob('data/sample_proteins/best_hits/*_best.txt')
seq_paths = 'data/reference_sequences/sporulation_sprot.fasta'

# parse the sequence descriptions to get a dict keyed by the first word
name_lookup = {}
for idx, gene in enumerate(SeqIO.parse(seq_paths, format='fasta')):
    spgene = gene.description.split(' ')
    if re.search(r"GN=([a-zA-Z0-9_]+) PE=", gene.description):
        species_id = spgene[0].split('|')[-1].split('_')[-1].lower()
        gene_name = re.findall(r"GN=([a-zA-Z0-9_]+) PE=", gene.description)[0] + '_' + species_id
        name_lookup[spgene[0]] = gene_name

# now we can go through the text files and do this thing
# we're going to make one big dataframe with columns [genome, gene, escore, bitscore]
row_list = []
for hp in hit_paths:
    genome = hp.split('/')[-1].replace('_best.txt', '')
    with open(hp, 'r') as hits:
        for line in hits.readlines():
            spline = line.strip().split('\t')
            row = {'genome': genome,
                   'gene': name_lookup[spline[1]],
                   'escore': float(spline[10]),
                   'bitscore': float(spline[11])}
            row_list.append(row)
    
gene_df = pd.DataFrame(row_list)
gene_df.to_csv('data/sample_proteins/best_hits/best_hits.csv', index=False)
