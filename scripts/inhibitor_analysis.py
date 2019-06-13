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
    df = blacklist(df)
    df = merge_multiple_inhibitors(df)
    sm = df.loc[:,'SMILES'].values.tolist()
    mw = list(map(lambda s: MolWt(fingerprinting.mol_to_smiles(s)), sm))
    df.loc[:,'mw'] = mw
    df.loc[:,'log_mw'] = np.log(np.array(mw))
    
    ''' Compare structures '''
#    bs = fingerprinting.fingerprint_from_smiles(sm, num_bits=64, max_path=9, method='topological')
#    bs = fingerprinting.fingerprint_from_smiles(sm, method='circular')
#    fingerprint_biplot(bs, fp_groups=df.loc[:,'EC'])
#    evaluate_global_separation(bs, dfm.loc[:,'EC'], plot=True)
#    print(evaluate_local_separation(bs, df.loc[:,'EC'], plot=True))
   
    ''' Grid search parameters for topological/circular '''
    grid_search_parameters(df, mode='circular', agg='median')
    
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
#    fingerprint_biplot(bs, fp_groups=df2.loc[:,'log_mw'])
#    evaluate_global_separation(bs, df2.loc[:,'EC'], plot=True)
    
def grid_search_parameters(df, mode, bit_range=[6,7,8,9,10,11,12,13],
                           radii=[1,2,3,4,5], path_range=[5,6,7,8,9,10],
                           agg='median', output='../data/params/'):
    ''' Computes average and global silhouette-like metrics for 
        different parameterizations of the topological or circular
        fingerprinting strategies '''
    sm = df.loc[:,'SMILES'].values.tolist()
    ec = df.loc[:,'EC'].values.tolist()
    output1 = output + mode + '_perf_local_' + agg +'s.csv'
    output2 = output + mode + '_perf_global_' + agg +'s.csv'
    eval_function = np.median if agg=='median' else np.mean
    
    if mode == 'topological':
        perf_local = np.zeros((len(bit_range), len(path_range)))
        perf_global = np.zeros((len(bit_range), len(path_range)))
        for i, bit_length in enumerate(bit_range):
            n_bits = 2 ** bit_length
            for j, max_path in enumerate(path_range):
                print('Testing n_bits =', n_bits, '| max_path =', max_path)
                bs = fingerprinting.fingerprint_from_smiles(sm, 
                    method='topological', max_path=max_path, num_bits=n_bits)
                local_vals = evaluate_local_separation(bs, ec, plot=False)
                global_vals = evaluate_global_separation(bs, ec, plot=False)
                perf_local[i,j] = eval_function(local_vals)
                perf_global[i,j] = eval_function(global_vals)
        df1 = pd.DataFrame(data=perf_local, index=2**np.array(bit_range), columns=path_range)
        df2 = pd.DataFrame(data=perf_global, index=2**np.array(bit_range), columns=path_range)
                
    elif mode == 'circular':
        perf_local = np.zeros((len(bit_range), len(radii)))
        perf_global = np.zeros((len(bit_range), len(radii)))
        for i, bit_length in enumerate(bit_range):
            n_bits = 2 ** bit_length
            for j, radius in enumerate(radii):
                print('Testing n_bits =', n_bits, '| radius =', radius)
                bs = fingerprinting.fingerprint_from_smiles(sm, 
                    method='circular', radius=radius, num_bits=n_bits)
                local_vals = evaluate_local_separation(bs, ec, plot=False)
                global_vals = evaluate_global_separation(bs, ec, plot=False)
                perf_local[i,j] = eval_function(local_vals)
                perf_global[i,j] = eval_function(global_vals)
        df1 = pd.DataFrame(data=perf_local, index=2**np.array(bit_range), columns=radii)
        df2 = pd.DataFrame(data=perf_global, index=2**np.array(bit_range), columns=radii)
    
    df1.to_csv(output1)
    df2.to_csv(output2)
    
    
def merge_multiple_inhibitors(df, bs=None):
    ''' Re-labels EC# of inhibitors that inhibit muliple ECs as "multiple" '''
    unique_smiles = set(); multiple_ec_smiles = set()
    for entry in df.index: # identify multiple EC inhibitors
        smiles = df.loc[entry,'SMILES']
        if smiles in unique_smiles: # multiple EC inhibitor
            multiple_ec_smiles.add(smiles)
        unique_smiles.add(smiles)
    counted_smiles = set(); duplicate_instances = []; unique_indices = []
    for i,entry in enumerate(df.index): # update EC field to say multiple
        smiles = df.loc[entry,'SMILES']
        if smiles in multiple_ec_smiles:
            df.loc[entry,'EC'] = 'multiple'
            if smiles in counted_smiles: # duplicate entry
                duplicate_instances.append(entry)
            else: # non-duplicate entry
                unique_indices.append(i)
            counted_smiles.add(smiles)
        else:
            unique_indices.append(i)
    dfm = df.drop(duplicate_instances)
    if bs is None:
        return dfm
    else:
        bsm = bs[np.array(unique_indices),:]
        return dfm,bsm
        
    
def blacklist(df, bs=None):
    ''' Remove metabolites if they meet any of the criteria:
        - Is inorganic (check by looking for carbon) 
        - Is variable (check by looking for * in SMILES string) '''
    remove = []; remove_indices = []
    for i, entry in enumerate(df.index):
        smiles = df.loc[entry,'SMILES']
        if '*' in smiles:
            remove.append(entry)
            remove_indices.append(i)
            print('Removed (variable):', smiles)
        else:
            carbon = fingerprinting.mol_to_smiles("C")
            mol = fingerprinting.mol_to_smiles(smiles)
            n_carbon = len(mol.GetSubstructMatches(carbon))
            if n_carbon == 0:
                remove.append(entry)
                remove_indices.append(i)
                print('Removed (inorganic):', smiles)
    n = df.shape[0]
    remove_indices = np.array(remove_indices)
    keep_indices = np.ones(n)
    np.put(keep_indices, remove_indices, 0)
    df_filt = df.drop(labels=remove)
    if bs is None:
        return df_filt
    else:
        bs_filt = bs[keep_indices.astype(bool),:]
        return df_filt, bs_filt
    
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
            max_diff = max(nearest_diff, nearest_same)
            if max_diff == 0:
                co = np.sign(nearest_diff - nearest_same)
            else:
                co = nearest_diff - nearest_same / max_diff
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


def fingerprint_biplot(fps, fp_groups=None, fp_labels=None):
    ''' 
    Visualize binary fingerprints in a biplot. Computes pairwise Jaccard 
    distances between each element, applies PCA, and projects the results 
    onto the top 2 PCs. Optionally, providing group labels will color-code
    the biplot based on group (such as for inhibitors, corresponding EC) 
    '''
    data = fingerprinting.fingerprint_jaccard_distances(fps, None, None)
    pca = sklearn.decomposition.PCA()
    top2 = pca.fit_transform(data)[:,:2]
    top2var = pca.explained_variance_ratio_[:2]
    X = top2[:,0]; Y = top2[:,1]
    
    if fp_groups is None:
        df = pd.DataFrame(data=[X,Y], index=['PC1','PC2']).T
        sp = sns.scatterplot(data=df, x='PC1', y='PC2')
    else:             
        df = pd.DataFrame(data=[X,Y,fp_groups], index=['PC1','PC2','group']).T
        sp = sns.scatterplot(data=df, x='PC1', y='PC2', hue='group')
        
    if not fp_labels is None:
        df.loc[:,'label'] = fp_labels
        xdiff = max(X) - min(X)
        for line in range(0,df.shape[0]):
            sp.text(df.loc[line,'PC1']+xdiff/100, df.loc[line,'PC2'], df.loc[line,'label'], 
                    horizontalalignment='left', size='small', color='black')
        
    plt.xlabel('PC1 (' + str(round(100*top2var[0],1)) + '%)')
    plt.ylabel('PC2 (' + str(round(100*top2var[1],1)) + '%)')
    plt.tight_layout()


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
    ax.axvline(x=np.median(silhouette_values), color="red", linestyle="--")
#        ax.text(silhouette_avg, y_lower + 1.5*y_shift, 'Median SC', ha='center')
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