#!/usr/bin/env python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier


# load the sample data
bitscores = pd.read_csv(
    'data/summary_data/sporulation_bitscores.csv')
genes = bitscores.drop('genome', axis='columns').columns
genomes = bitscores['genome']

# load the classes
class_df = pd.read_csv('data/sample_classes.csv', index_col=0)
bitscores = pd.merge(bitscores, class_df, on='genome')
labels = bitscores['class'].values

# pull out and split the data
mat = bitscores.drop(['genome', 'class'], axis='columns').values
pa_mat = (mat > 0).astype(int)


# We're going to populate a vector of gene importances across all folds
gene_importances = np.zeros(genes.shape[0])

# get fold labels for our k-fold model fits
k_folds = 10
fold_labels = np.zeros(genomes.shape[0])
repeat_nums = np.repeat(range(k_folds), int(fold_labels.shape[0] / k_folds))
fold_labels[:repeat_nums.shape[0]] = repeat_nums
rng = np.random.default_rng(seed=2046)
fold_labels = rng.permutation(fold_labels)

for ki in np.arange(k_folds):
    # do the train test split
    test_idx = fold_labels == ki
    train_idx = fold_labels != ki
    X_train = pa_mat[train_idx, :]
    y_train = labels[train_idx]
    X_test = pa_mat[test_idx, :]
    y_test = labels[test_idx]

    # fit the random forest classifier
    rfc = RandomForestClassifier(max_depth=1, n_estimators=1000)
    rfc.fit(X_train, y_train)
    y_pred = rfc.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    print(f"Fold 0 Cluster Recovery: {acc*100:.2f} %")

    # add gene importances
    gene_importances += rfc.feature_importances_

    # if its the first fold we'll plot it
    if ki == 0:
        # We're going to try visualizing in PCA but classifying in real space
        pca = PCA()
        X_train_pca = pca.fit_transform(X_train)
        X_test_pca = pca.transform(X_test)

        # need grid of predictions
        resolution = 50
        xmin = np.min(X_train_pca[:, 0]) - 1
        xmax = np.max(X_train_pca[:, 0]) + 1
        ymin = np.min(X_train_pca[:, 1]) - 1
        ymax = np.max(X_train_pca[:, 1]) + 1
        x_range = np.linspace(xmin, xmax, resolution)
        y_range = np.linspace(ymin, ymax, resolution)
        xx, yy = np.meshgrid(x_range, y_range)

        # make the predictions
        preds = np.zeros(xx.shape)
        pca_feature_vec = np.zeros(X_train_pca.shape[1])
        for xi, x in enumerate(x_range):
            for yi, y in enumerate(y_range):
                pca_feature_vec[0], pca_feature_vec[1] = xx[xi, yi], yy[xi, yi]
                feature_vec = pca.inverse_transform(
                    pca_feature_vec.reshape(1, -1))
                preds[xi, yi] = rfc.predict(feature_vec)

        # lets look at the classification
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.contourf(xx, yy, preds, levels=2, cmap='RdBu', alpha=0.3)
        ax.scatter(X_test_pca[:, 0], X_test_pca[:, 1],
                   c='k', marker='x', alpha=0.3, s=30)

        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f} %)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f} %)")
        ax.set_aspect('equal')
        ax.set_title(f"RFC Cluster Recovery: {acc * 100:.3f}%")
        plt.savefig('plots/raw/rfc.png')

out_df = pd.DataFrame({'gene': genes, 'importance': gene_importances})
out_df.to_csv('data/summary_data/rfc_genes.csv')
