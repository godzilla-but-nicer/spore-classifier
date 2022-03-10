#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import time

# load the gene list data

mode = sys.argv[1]

if mode == "sample":
    prot_dir = "data/sample_proteins"
elif mode == "ww":
    prot_dir = "data/labelled/ww_proteins"
elif mode == "browne":
    prot_dir = "data/labelled/browne_proteins"

genes_df = pd.read_csv(prot_dir + '/best_hits/best_hits.csv')

# pivot it to a pres/abs matrix with
genes_df['bitscore'] = genes_df['bitscore'].astype(float)
mat = (genes_df
      .drop_duplicates(['genome', 'gene']) # this is hacky, we should always keep higher bitscore
      .drop('escore', axis='columns')
      .pivot(index='genome', columns='gene', values='bitscore')
      .fillna(0.0))

# save the presence absence matrix
if mode == "sample":
    mat.to_csv('data/summary_data/sporulation_bitscores.csv')
elif mode == "ww":
    mat.to_csv('data/summary_data/ww_sporulation_bitscores.csv')
elif mode == "browne":
    mat.to_csv('data/summary_data/browne_sporulation_bitscores.csv')


plt.imshow(mat.values, interpolation='none', cmap='Greys')
plt.savefig('plots/raw/bitscore_mat.png')