#!/usr/bin/env python

import os
import pandas as pd
import sys

# load the data from the browne et al paper
browne = pd.read_csv('data/labelled/browne_et_al_labels.csv', skiprows=2, index_col=0)

# first of all we need to make some directories if they dont exist
browne_dir = "data/labelled/browne_proteins"
if not os.path.isdir(browne_dir):
    sys.mkdir(browne_dir)

for num in browne['accession_number']:
    # if we're given an accession number we have to do some stuff
    if num != 'cultured in this study':
        os.system(f"./datasets download genome accession {num} " 
                 + "--exclude-gff3 --exclude-rna "
                 + "--exclude-seq --exclude-genomic-cds")
        os.system("unzip -qo ncbi_dataset.zip")
        os.system(f"mv ncbi_dataset/data/{num}/protein.faa {browne_dir}/{num}.faa")
        os.system("rm -r ncbi_dataset/")


