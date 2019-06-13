#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 19:09:34 2019

@author: jhyun95
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.decomposition

from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.Chem.Fingerprints.FingerprintMols import FingerprintMol, FingerprintsFromMols # topological
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect # circular

def fingerprint_from_smiles(smiles, method='topological',
        min_path=1, max_path=7, radius=2, num_bits=2048):
    ''' 
    Generates a chemical fingerprint of a compound from a SMILES string
    
    Parameters 
    ----------
    smiles : str or list
        SMILES string(s) to convert to binary fingerprint(s)
    method : str
        "topological" or "daylight" for topological fingerprinting, "circular" 
        or "morgan" for circular fingerprinting (default "topological)
    min_path : int
        Minimum path length for topological fingerprinting (default 1)
    max_path : int
        Maximum path length for topological fingerprinting (default 7)
    radius : int
        Fingerprint radius for circular fingerprinting (default 2)
    num_bits : int
        Length of fingerprint bit vector, for both methods (default 2048)
        
    Returns
    -------
    out : ndarray
        Numpy integer array with fingerprint bit vector
    '''
    
    def fingerprint_to_array(fp):
        ''' Convert rdkit ExplicitBitVect to numpy boolean array
            Note: Was unable to make use of fp.ToBinary() could not reliably
            decode the returned byte string into the correct bit string. See:
            https://sourceforge.net/p/rdkit/mailman/message/24426410/ '''
        bs = fp.ToBitString() 
        bs_list = list(map(int,bs))
        return np.array(bs_list, dtype=int)
        
    if type(smiles) == str: # single entry
        mol = mol_to_smiles(smiles)
        if method=='topological' or method=='daylight':
            # https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RDKFingerprint
            fp = FingerprintMol(mol, minPath=min_path, maxPath=max_path, 
                fpSize=num_bits, bitsPerHash=2, useHs=True, tgtDensity=0, minSize=128)
        elif method=='circular' or method=='morgan':
            fp = GetMorganFingerprintAsBitVect(mol, radius, nBits=num_bits)
        return fingerprint_to_array(fp)
    
    elif type(smiles) == list: # multiple entries
        mols = map(lambda i: (i, mol_to_smiles(smiles[i])), range(len(smiles)))
        if method=='topological' or method=='daylight':
            fps = FingerprintsFromMols(list(mols), minPath=min_path, maxPath=max_path, 
                fpSize=num_bits, bitsPerHash=2, useHs=True, tgtDensity=0, minSize=128,
                reportFreq=-1)
            fps = map(lambda x: x[1], fps)
        elif method=='circular' or method=='morgan':
            fps = map(lambda mol: GetMorganFingerprintAsBitVect(mol[1], 
                radius, nBits=num_bits), mols)
        bss = list(map(fingerprint_to_array, fps))
        return np.array(bss)
    

def fingerprint_jaccard_distances(fps, fp_labels=None, visualize=None):
    ''' 
    Computes pairwise Jaccard distances between binary fingerprints,
    and optionally creates a heatmap or clustermap to visualize distances
    
    Parameters
    ----------
    fps : ndarray
        2D integer or boolean Numpy array representing chemical x fingerprint
    fp_labels : list
        List of strings with chemical labels for visualization (default None)
    visualize : str
        Type of visualization to generate (either 'heatmap' or 'clustermap'), 
        or other value to skip visualization (default None).
        
    Returns
    -------
    distances : ndarray
        Symmetric Numpy float array with pairwise distances
    '''
    
    n = fps.shape[0]; distances = np.zeros((n,n))
    jaccard = lambda x,y: np.logical_and(x,y).sum() / float(np.logical_or(x,y).sum())
    for i in range(n):
        for j in range(i):
            dist = 1.0 - jaccard(fps[i,:], fps[j,:])
            if pd.isnull(dist): # NaN, if both vectors are all 0s
                dist = 1.0
            distances[i,j] = dist
            distances[j,i] = dist
    
    labels = range(n) if fp_labels is None else fp_labels
    df_dists = pd.DataFrame(data=distances, index=labels, columns=labels)
    if visualize == 'clustermap': # optionally make a clustermap
        sns.clustermap(df_dists)
    elif visualize == 'heatmap': # or a heatmap
        sns.heatmap(df_dists, vmin=0.0, vmax=1.0)
    return distances
    

def mol_to_smiles(smiles_str):
    ''' Wrapper for MolFromSmiles to avoid getting stuck on abnormal strings '''
    mol = MolFromSmiles(smiles_str)
    if mol is None: # failed to run with sanitization
        print('SMILES', smiles_str, 'is abnormal, skipping sanitization')
        mol = MolFromSmiles(smiles_str, sanitize=False)
        mol.UpdatePropertyCache(strict=False)
    return mol