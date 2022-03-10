#!/bin/usr/env python

import pandas as pd
from sklearn.model_selection import train_test_split

bitscores = pd.read_csv('data/sample_classes.csv', index_col=0)
genomes_select, genomes_verify, \
    labels_select, labels_verify = train_test_split(
        bitscores['genome'].values, bitscores['class'].values)
