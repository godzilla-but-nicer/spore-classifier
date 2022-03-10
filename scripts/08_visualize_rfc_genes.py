#!/bin/usr/env python

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

rfc_genes = pd.read_csv('data/summary_data/rfc_genes.csv', index_col=0)
nonzero = rfc_genes[rfc_genes['importance'] > 0]
rfc_nonzero = nonzero.sort_values('importance', ascending=False)
rfc_nonzero['importance'] /= 10
rfc_nonzero.to_csv('data/summary_data/rfc_nonzero.csv')

# get the cumulative importance
cumulative_importance = np.cumsum(rfc_nonzero['importance'])


# make the plot
fig, ax = plt.subplots()
#ax.scatter(rfc_nonzero['gene'], rfc_nonzero['importance']*10)
ax.plot(rfc_nonzero['gene'], cumulative_importance, marker='o', markersize=2)
ax.set_xticklabels([], rotation=90, size=5)
ax.set_xticks([])
ax.set_xlabel('Ranked Genes')
ax.set_ylabel('RFC Importance')
plt.show()
