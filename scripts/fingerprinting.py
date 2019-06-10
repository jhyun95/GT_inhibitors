#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 19:09:34 2019

@author: jhyun95
"""

import numpy as np
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
        mol = MolFromSmiles(smiles)
        if method=='topological' or method=='daylight':
            # https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RDKFingerprint
            fp = FingerprintMol(mol, minPath=min_path, maxPath=max_path, 
                fpSize=num_bits, bitsPerHash=2, useHs=True, tgtDensity=0, minSize=128)
        elif method=='circular' or method=='morgan':
            fp = GetMorganFingerprintAsBitVect(mol, radius, nBits=num_bits)
        return fingerprint_to_array(fp)
    
    elif type(smiles) == list: # multiple entries
        mols = map(lambda i: (i, MolFromSmiles(smiles[i])), range(len(smiles)))
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
        