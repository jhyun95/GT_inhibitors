# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:37:48 2019

@author: jhyun_000
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import sklearn.metrics, sklearn.preprocessing

import fingerprinting
from rdkit.Chem.Descriptors import MolWt

def main():
    ''' Load structures and compute MWs '''
    df = pd.read_csv('../data/inhibitor_merged.tsv', delimiter='\t')
    sm = df.loc[:,'SMILES'].values.tolist()
    mw = list(map(lambda s: MolWt(fingerprinting.mol_to_smiles(s)), sm))
    df.loc[:,'mw'] = mw
    df.loc[:,'log_mw'] = np.log(np.array(mw))
    
    ''' Compare structures '''
#    bs = fingerprinting.fingerprint_from_smiles(sm, method='topological')
    bs = fingerprinting.fingerprint_from_smiles(sm, method='circular')
#    fingerprinting.fingerprint_biplot(bs, fp_groups=df.loc[:,'log_mw'])
    evaluate_global_separation(bs, df.loc[:,'EC'], plot=True)
#    print(evaluate_local_separation(bs, df.loc[:,'EC'], plot=True))
    
    ''' Match SNG output existing tables '''
#    df_sng = pd.read_csv('../data/scaffold_fingerprints_inhibitors.csv', index_col=0).T
#    df_sng.index = df_sng.index.astype(int)
#    sng_dim = df_sng.shape[1] # size of SNG fingerprint
#    df2 = df.copy()
#    for col in df_sng.columns:
#        df2[col] = np.nan
#    for i in range(df2.shape[0]): # this loop is really slow 
#        gid = df2.loc[i,'group_id']
#        if gid in df_sng.index:
#            df2.iloc[i,-sng_dim:] = df_sng.loc[gid,:].values
#    df2 = df2.dropna(how='any')
#    bs = df2.iloc[:,-sng_dim:].values
#    fingerprinting.fingerprint_biplot(bs, fp_groups=df2.loc[:,'log_mw'])
#    evaluate_global_separation(bs, df2.loc[:,'EC'], plot=True)
    
def evaluate_local_separation(fps, ecs, plot=False):
    ''' Evaluate local separation by computing the ratio of the distance 
        to the nearest point in the same EC# vs. nearest point in a
        different EC#. '''
    distances = fingerprinting.fingerprint_jaccard_distances(fps)
    le = sklearn.preprocessing.LabelEncoder(); le.fit(ecs)
    ec_labels = le.transform(ecs)
    
    cohesiveness = np.zeros(fps.shape[0])
    for i in range(fps.shape[0]):
        current_ec = ec_labels[i]
        same_ec_entries = (ec_labels==current_ec)
        diff_ec_entries = np.logical_not(same_ec_entries)
        same_ec_entries[i] = False # ignore self
        if np.sum(same_ec_entries.astype(int)) == 0: # singleton cluster
            cohesiveness[i] = -1.0
        else:
            dist_to_same = distances[i,same_ec_entries]
            dist_to_diff = distances[i,diff_ec_entries]
            nearest_same = np.min(dist_to_same)
            nearest_diff = np.min(dist_to_diff)
            co = nearest_diff - nearest_same / max(nearest_diff, nearest_same)
            cohesiveness[i] = co
    
    if plot: # generate silhouette plot
        ax = plot_silhouette(cohesiveness, ecs, metric_label='Local cohesiveness')
    return cohesiveness
    
def evaluate_global_separation(fps, ecs, plot=False):
    ''' Evaluate global separating by computing silhouette coefficients, 
        treating EC#s as clusters '''
    distances = fingerprinting.fingerprint_jaccard_distances(fps)
    silhouette_values = sklearn.metrics.silhouette_samples(distances, ecs, metric='precomputed')
    silhouette_values[np.isnan(silhouette_values)] = -1.0 # replace nan with worst score
    if plot: # generate silhouette plot
        ax = plot_silhouette(silhouette_values, ecs)
    return silhouette_values


def plot_silhouette(silhouette_values, labels, group_label='EC Number', 
                    metric_label="Silhouette Coefficient (SC)"):
    # https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html
    le = sklearn.preprocessing.LabelEncoder(); le.fit(labels)
    ec_labels = le.transform(labels)
    unique_ecs = list(le.classes_)
    
    fig, ax = plt.subplots(1,1)
    y_shift = 10; y_lower = y_shift
    for i in range(len(unique_ecs)):
        ith_cluster_silhouette_values = silhouette_values[ec_labels==i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        color = cm.nipy_spectral(float(i) / len(unique_ecs))
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, unique_ecs[i])
        y_lower = y_upper + y_shift 
    
    ax.set_xlabel(metric_label)
    ax.set_ylabel(group_label)
    ax.axvline(x=np.mean(silhouette_values), color="red", linestyle="--")
#        ax.text(silhouette_avg, y_lower + 1.5*y_shift, 'Average SC', ha='center')
    ax.set_yticks([])
    return ax


def merge_inhibitor_data(output='../data/inhibitor_merged.tsv'):
    ''' Merge tables with inhibitor SMILES and common names'''
    df1 = pd.read_csv('../data/inhibitor_smiles.tsv', delimiter='\t')
    df2 = pd.read_csv('../data/inhibitor_names.csv')
    df1['name'] = None
    n = df1.shape[0]; m = df2.shape[0]
    
    ''' Get all molecule names for (ec, group_id) pairs '''
    entry_to_name = {}
    for i in range(m): # get all (ec,group_id) name labels
        ec = df2.loc[i,'EC']
        gid = df2.loc[i,'group_id']
        entry = (ec,gid)
        name = df2.loc[i,'name']
        if not entry in entry_to_name: 
            entry_to_name[entry] = name
        elif entry_to_name[entry] == 'more': # replace duplicates that just say "more"
            entry_to_name[entry] = name
            
    ''' Get all SMILES for (ec, group_id) pairs, map SMILES:name '''
    smiles_to_name = {}
    for i in range(n): # merge results into smiles table, attempt to resolve dupes
        ec = df1.loc[i,'EC']; 
        gid = df1.loc[i,'group_id']
        smiles = df1.loc[i,'SMILES']
        if (ec,gid) in entry_to_name:
            name = entry_to_name[(ec,gid)]
            if name != 'more': # not a duplicate
                smiles_to_name[smiles] = name
    
    ''' Fill in names, prioritizing SMILES:name mapping to avoid duplicates ''' 
    df1.loc[:,'name'] = None
    for i in range(n): # merge results into smiles table, attempt to resolve dupes
        ec = df1.loc[i,'EC']; 
        gid = df1.loc[i,'group_id']
        smiles = df1.loc[i,'SMILES']
        name = None
        if smiles in smiles_to_name:
            name = smiles_to_name[smiles]
        elif (ec,gid) in entry_to_name:
            name = entry_to_name[ec,gid]
        df1.loc[i,'name'] = name
    
    if not output is None:
        df1.to_csv(output, sep='\t')
    
    return df1
            
    
if __name__ == '__main__':
    main()