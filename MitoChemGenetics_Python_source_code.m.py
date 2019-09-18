''' Compute the p-values and FDRs for a list of Z-scores or deltaZscores'''
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import *
from scipy.stats import norm
from statsmodels.stats import multitest

''' For a list of Z-scores'''
Z = pd.read_csv('../data/z_score.csv', index_col=0)
drugs = Z.columns.tolist()
dfs = {}
for drug in drugs:
    z = Z[drug]
    p = norm.cdf(-z.abs()) + (1 - norm.cdf(z.abs()))
    f = multitest.multipletests(p, method='fdr_bh')[1]
    p = -np.log10(p)
    df = pd.DataFrame({'Z': z, ' -log10(p-val)': p, 'FDR': f}).sort_values('FDR', ascending=True)
    dfs[drug] = df
    df.to_csv('../figures/new_Z_%s.csv' % drug.replace('/', '_'))

''' For a list of deltaZscores'''
Z = pd.read_csv('../data/delta_z_score.csv', index_col=0)
drugs = Z.columns.tolist()
dfs = {}
for drug in drugs:
    z = Z[drug]
    p = norm.cdf(-z.abs()) + (1 - norm.cdf(z.abs()))
    f = multitest.multipletests(p, method='fdr_bh')[1]
    p = -np.log10(p)
    df = pd.DataFrame({'Z': z, ' -log10(p-val)': p, 'FDR': f}).sort_values('FDR', ascending=True)
    dfs[drug] = df
    df.to_csv('../figures/new_delta_Z_%s.csv' % drug.replace('/', '_'))
