#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA 
from sklearn.cluster import KMeans
from umap import UMAP

# load the sample data
bitscores = pd.read_csv('data/summary_data/sporulation_bitscores.csv')
genes = bitscores.drop('genome', axis='columns').columns
genomes = bitscores['genome']

# load the weller and wu data



# Do the dimensionality reduction
mat = bitscores.drop('genome', axis='columns').values
pca = PCA(n_components=10)
umap = UMAP()

# let's first get a distribution of row sums for the prescence
pa_mat = (mat > 0).astype(int)
row_sums = np.sum(pa_mat, axis=1)
plt.figure()
plt.hist(row_sums, bins=30)
plt.savefig('plots/raw/genecount_hist.png')


# bitscore pca
mat_pca = pca.fit_transform(mat)

plt.figure()
plt.scatter(mat_pca[:,0], mat_pca[:,1], s=10)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f} %)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f} %)")
plt.savefig('plots/raw/pca_bitscores.png')


# presence absence PCA

pa_mat_pca = pca.fit_transform(pa_mat)

plt.figure()
plt.scatter(pa_mat_pca[:,0], pa_mat_pca[:,1], s=10)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f} %)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f} %)")
plt.savefig('plots/raw/pca_prescence_absence.png')

# fit umap and visualize
mat_umap = umap.fit_transform(mat)

plt.figure()
plt.scatter(mat_umap[:,0], mat_umap[:,1], s=10)
plt.xlabel(f"UMAP 1")
plt.ylabel(f"UMAP 2")
plt.savefig('plots/raw/umap_bitscore.png')

# umap on prescence absence
pa_mat_umap = umap.fit_transform(mat)

plt.figure()
plt.scatter(pa_mat_umap[:,0], pa_mat_umap[:,1], s=10)
plt.xlabel(f"UMAP 1")
plt.ylabel(f"UMAP 2")
plt.savefig('plots/raw/umap_prescence_absence.png')