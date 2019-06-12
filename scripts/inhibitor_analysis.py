# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:37:48 2019

@author: jhyun_000
"""

import pandas as pd
import numpy as np
import seaborn as sns
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
    bs= fingerprinting.fingerprint_from_smiles(sm, method='circular')
    fingerprinting.fingerprint_biplot(bs, fp_groups=df.loc[:,'log_mw'])
    
#    df_sng = pd.read_csv('../data/scaffold_fingerprints_inhibitors.csv', index_col=0)
    
    
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