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

# Do the dimensionality reduction on the pres/abs data
mat = bitscores.drop('genome', axis='columns').values
pa_mat = (mat > 0).astype(int)
pca = PCA(n_components=2)
pa_mat_pca = pca.fit_transform(pa_mat)

# do the clustering
km = KMeans(n_clusters=2)
km.fit(pa_mat_pca)

# predict over a grid for the plot
# grid setup
resolution = 200
xmin, xmax = (np.min(pa_mat_pca[:, 0]), np.max(pa_mat_pca[:, 0]))
xspace = np.linspace(xmin-1, xmax+1, resolution)
ymin, ymax = (np.min(pa_mat_pca[:, 1]), np.max(pa_mat_pca[:, 1]))
yspace = np.linspace(ymin-1, ymax+1, resolution)
xx, yy = np.meshgrid(xspace, yspace)

# predictions
preds = np.zeros(xx.shape)
for xi, _ in enumerate(xspace):
    for yi, _ in enumerate(yspace):
        preds[xi, yi] = km.predict(
            np.array((xx[xi, yi], yy[xi, yi])).reshape(1, -1))

# plot the prediction regions
fig, ax = plt.subplots()
ax.contourf(xx, yy, preds, cmap='Accent')
ax.scatter(pa_mat_pca[:, 0], pa_mat_pca[:, 1], c='k', marker='x')
ax.set_aspect('equal')
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f} %)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f} %)")
plt.savefig('plots/raw/2_clusters.png')

# make the cluster predictions for the two classes
# and save a df with this info
classes = km.predict(pa_mat_pca)
class_df = pd.DataFrame({'genome': genomes, 'class': classes})
class_df.to_csv('data/sample_classes.csv')


# just for funsies
km = KMeans(n_clusters=3)
km.fit(pa_mat_pca)

# predict over a grid for the plot
# predictions
preds = np.zeros(xx.shape)
for xi, _ in enumerate(xspace):
    for yi, _ in enumerate(yspace):
        preds[xi, yi] = km.predict(
            np.array((xx[xi, yi], yy[xi, yi])).reshape(1, -1))

# plot the prediction regions
fig, ax = plt.subplots()
ax.contourf(xx, yy, preds, cmap='Accent')
ax.scatter(pa_mat_pca[:, 0], pa_mat_pca[:, 1], c='k', marker='x')
ax.set_aspect('equal')
ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f} %)")
ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f} %)")
plt.savefig('plots/raw/3_clusters.png')
#ww_pca = pca.transform(ww_pa_mat)
#browne_pca = pca.transform(browne_pa_mat)
#label_convert = ['N', 'Y']
#ww_colors = ['red', 'blue']
# for i in range(2):
#    #plot_idx = np.where(cluster_ids == i)
#    ww_idx = np.where(ww_labels == label_convert[i])
#    browne_idx = np.where(browne_labels == i)
#    # ww points
#    plt.scatter(ww_pca[ww_idx, 0], ww_pca[ww_idx, 1],
#                label=f'Weller and Wu: {bool(i)}', facecolor='none',
#                edgecolor=ww_colors[i], s=12)
#    plt.scatter(browne_pca[browne_idx, 0], browne_pca[browne_idx, 1],
#                label=f'Browne et al.: {bool(i)}', facecolor='none',
#                edgecolor=ww_colors[i], marker='s', s=12)
#
