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
ww_bitscores = pd.merge(ww_bitscores, ww_labels_df, on='genome')

# we need to make sure that both dfs have the same columns
ww_genes = ww_bitscores.drop(
    ['genome', 'spore_forming'], axis='columns').columns
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

# Ok we have to do the same shit but for the browne data
browne_bitscores = pd.read_csv(
    'data/summary_data/browne_sporulation_bitscores.csv')
browne_labels = pd.read_csv(
    'data/labelled/browne_et_al_labels.csv', skiprows=2, index_col=0)
browne_labels['genome'] = browne_labels['accession_number']
browne_bitscores = pd.merge(browne_bitscores, browne_labels[[
                            'genome', 'category']], on='genome', how='left')

# make sure they all have the same columns
browne_genes = browne_bitscores.drop(
    ['genome', 'category'], axis='columns').columns
missing_cols = np.array([g for g in genes if g not in browne_genes])
for col in missing_cols:
    browne_bitscores[col] = 0

# here we're dropping extra columns that dont show up in our unlabelled data
keep_cols = list(bitscores.columns)
keep_cols.append('category')
keep_cols = np.array(keep_cols)
browne_bitscores = browne_bitscores[keep_cols]

# get the objects we will actually use
browne_mat = browne_bitscores.drop(
    ['genome', 'category'], axis='columns').values
browne_labels_string = browne_bitscores['category'].values
browne_labels = np.array(
    [0 if "non-spore" in lab else 1 for lab in browne_labels_string])
browne_pa_mat = (browne_mat > 0).astype(int)

# Do the dimensionality reduction on the bitscore data
mat = bitscores.drop('genome', axis='columns').values
pca = PCA(n_components=10)

# bitscore pca
mat_pca = pca.fit_transform(mat)

plt.figure()
plt.scatter(mat_pca[:, 0], mat_pca[:, 1], s=10)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f} %)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f} %)")
plt.savefig('plots/raw/pca_bitscores.png')


# presence absence PCA
pa_mat = (mat > 0).astype(int)
pa_mat_pca = pca.fit_transform(pa_mat)

plt.figure()
plt.scatter(pa_mat_pca[:, 0], pa_mat_pca[:, 1], s=10)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f} %)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f} %)")
plt.savefig('plots/raw/pca_prescence_absence.png')

# fit umap and visualize
umap = UMAP(n_neighbors=50)
mat_umap = umap.fit_transform(mat)

plt.figure()
plt.scatter(mat_umap[:, 0], mat_umap[:, 1], s=10)
plt.xlabel(f"UMAP 1")
plt.ylabel(f"UMAP 2")
plt.savefig('plots/raw/umap_bitscore.png')

# umap on prescence absence
pa_mat_umap = umap.fit_transform(pa_mat)

plt.figure()
plt.scatter(pa_mat_umap[:, 0], pa_mat_umap[:, 1], s=10)
plt.xlabel(f"UMAP 1")
plt.ylabel(f"UMAP 2")
plt.savefig('plots/raw/umap_prescence_absence.png')

# pca the weller and wu data
ww_pca = pca.transform(ww_pa_mat)
browne_pca = pca.transform(browne_pa_mat)

fig, ax = plt.subplots(figsize=(8, 8))
plt.scatter(pa_mat_pca[:, 0], pa_mat_pca[:, 1],
            s=10, c='k', marker='x', alpha=0.2)
# ww_labels are in a char format
label_convert = ['N', 'Y']
ww_colors = ['red', 'blue']
for i in range(2):
    #plot_idx = np.where(cluster_ids == i)
    ww_idx = np.where(ww_labels == label_convert[i])
    browne_idx = np.where(browne_labels == i)
    # ww points
    ax.scatter(ww_pca[ww_idx, 0], ww_pca[ww_idx, 1],
               label=f'Weller and Wu: {bool(i)}', facecolor='none',
               edgecolor=ww_colors[i], s=16)
    ax.scatter(browne_pca[browne_idx, 0], browne_pca[browne_idx, 1],
               label=f'Browne et al.: {bool(i)}', facecolor='none',
               edgecolor=ww_colors[i], marker='s', s=16)


ax.set_aspect('equal')


ax.legend()
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f} %)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f} %)")
plt.savefig('plots/raw/labelled.png')
