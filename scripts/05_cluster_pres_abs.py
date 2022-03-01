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
ww_bitscores = pd.read_csv('data/summary_data/ww_sporulation_bitscores.csv')
ww_labels_df = pd.read_csv(
    'data/labelled/ww_genomes_labelled.csv').drop(['umap_0', 'umap_1'], axis='columns')
ww_bitscores = pd.merge(ww_bitscores, ww_labels_df, on = 'genome')

# we need to make sure that both dfs have the same columns
ww_genes = ww_bitscores.drop(['genome', 'spore_forming'], axis = 'columns').columns
missing_cols = np.array([g for g in genes if g not in ww_genes])
for col in missing_cols:
    ww_bitscores[col] = 0

keep_cols = list(bitscores.columns)
keep_cols.append('spore_forming')
keep_cols = np.array(keep_cols)
ww_bitscores = ww_bitscores[keep_cols]

# get a bunch of useful vectors from ww data
ww_genomes = ww_bitscores['genome']
ww_mat = ww_bitscores.drop(['genome', 'spore_forming'], axis='columns').values
ww_labels = ww_bitscores['spore_forming'].values
ww_pa_mat = (ww_mat > 0).astype(int)

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
# mat_umap = umap.fit_transform(mat)
# 
# plt.figure()
# plt.scatter(mat_umap[:,0], mat_umap[:,1], s=10)
# plt.xlabel(f"UMAP 1")
# plt.ylabel(f"UMAP 2")
# plt.savefig('plots/raw/umap_bitscore.png')
# 
# # umap on prescence absence
# pa_mat_umap = umap.fit_transform(mat)
# 
# plt.figure()
# plt.scatter(pa_mat_umap[:,0], pa_mat_umap[:,1], s=10)
# plt.xlabel(f"UMAP 1")
# plt.ylabel(f"UMAP 2")
# plt.savefig('plots/raw/umap_prescence_absence.png')

# pca the weller and wu data
ww_pca = pca.transform(ww_pa_mat)

# setup kmeans and run on the PCA data
km = KMeans(n_clusters=2)
cluster_ids = km.fit_predict(pa_mat_pca[:, :2])

# get a grid of predictions
buffer = 1
resolution = 1000
x = np.linspace(np.min(pa_mat_pca[:,0]) - buffer, np.max(pa_mat_pca[:, 0]) + buffer, resolution)
y = np.linspace(np.min(pa_mat_pca[:,1]) - buffer, np.max(pa_mat_pca[:, 1]) + buffer, resolution)
xx, yy = np.meshgrid(x, y)
predictions = np.zeros(xx.shape)
for xi, _ in enumerate(x):
    for yi, _ in enumerate(y):
        predictions[xi, yi] = km.predict(np.array([x[xi], y[yi]]).reshape(1, -1))



plt.figure()
plt.contourf(xx, yy, predictions.T, alpha = 0.2)
# sample points
plt.scatter(pa_mat_pca[:, 0], pa_mat_pca[:, 1], s=10, c='k', marker = 'x', alpha=0.4)
# ww_labels are in a char format
label_convert = ['N', 'Y']
ww_colors = ['red', 'blue']
for i in range(2):
    plot_idx = np.where(cluster_ids == i)
    ww_idx = np.where(ww_labels == label_convert[i])
    # ww points
    plt.scatter(ww_pca[ww_idx, 0], ww_pca[ww_idx, 1], c=ww_colors[i], label=f'Spore-forming: {bool(i)}')

plt.legend()
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f} %)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f} %)")
plt.savefig('plots/raw/clustered.png')
