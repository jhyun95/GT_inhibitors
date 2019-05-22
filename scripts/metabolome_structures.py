#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 14:47:27 2019

@author: jhyun95
"""

import os, urllib
import cobra
import pandas as pd
import rdkit.Chem

MODEL_FILE = '../data/iCHOv1.json'
METABOLITES_FILE = '../data/iCHOv1.chemicals.tsv'
METABOLITE_TAG_FILE = '../data/iCHOv1.chemical_tags.csv'
METABOLITE_MOL_DIR = '../raw_data/CHO_structures/'
METABOLITE_SMILES_FILE = '../data/iCHOv1_smiles.tsv'
INHIBITOR_DIR = '../raw_data/inhibitors/'
INHIBITOR_SMILES_FILE = '../data/inhibitor_smiles.tsv'

NAME_FIXES = { # for BiGG IDs that are on MetaNetX but not iCHOv1?
    'dhnpt': '7,8-dihydroneopterin',
    'sdhlam': 'S-succinyldihydrolipoamide',
    'ficytc': 'ferricytochrome',
    'pepd': 'peptide',
    'm2macchitppdol': '(alpha-D-mannosyl)2-beta-D-mannosyldiacetylchitobiosyldiphosphodolichol',
    'm1macchitppdol': 'alpha-D-mannosyl-beta-D-mannosyl-diacylchitobiosyldiphosphodolichol',
    'm3macchitppdol': '(alpha-D-mannosyl)3-beta-D-mannosyl-diacetylchitodiphosphodolichol',
    'ptd2meeta_SC': 'phosphatidyldi-N-methylethanolamine_tomerge'}

STRUCTURE_BLACKLIST = {'pepd'} # has X in structure, ambiguous peptide

db_priorities = ['chebi', 'hmdb', 'kegg', 'lipidmaps'] # structures readily available
db_priorities += ['seed', 'metacyc', 'sabiork', 'reactome', 'mnxm', 'bigg'] # structures less available

def main():
    ''' Downloading MOL files for CHO metabolites '''
    #df = process_metanetx()
    #download_mol_files()
    
    ''' Converting CHO metabolite structures to SMILES '''
#    mol_files = os.listdir(METABOLITE_MOL_DIR)
#    mol_paths = list(map(lambda x: METABOLITE_MOL_DIR + x, mol_files))
#    mol_to_smiles(mol_paths, smiles_path=METABOLITE_SMILES_FILE)
    
    ''' Converting GT inhibitor structures to SMILES '''
    with open(INHIBITOR_SMILES_FILE, 'w+') as f:
        f.write('EC\tgroup_id\tSMILES\n')
        for EC_num in os.listdir(INHIBITOR_DIR): # for each EC number
            EC_path = INHIBITOR_DIR + EC_num + '/'
            if os.path.isdir(EC_path):
                mol_paths = []
                for group_id in os.listdir(EC_path): # for each compound/group id
                    mol_path = EC_path + group_id
                    if not '.DS_Store' in mol_path:
                        mol_paths.append(mol_path)
                smiles = mol_to_smiles(mol_paths)
                for mol_name in smiles:
                    smiles_str = smiles[mol_name]
                    f.write(EC_num + '\t' + mol_name + '\t' + smiles_str + '\n')
    
    
def mol_to_smiles(mol_path, smiles_path=None):
    ''' Converts a MOL file(s) to SMILES format. Returns a dictionary
        with keys as the molecule name (filepath with .mol removed), and
        values as the SMILES strings. Optionally output to file (tsv) '''
    if type(mol_path) == str: # single file
        mol_files = [mol_path]
    else: # multiple files
        mol_files = mol_path
        
    ''' Convert each MOL file to SMILES strings '''
    smiles = {}
    for mol_file in mol_files:
        print('Converting', mol_file)
        mol_name = mol_file.split('/')[-1][:-4]
        if not mol_name in STRUCTURE_BLACKLIST:
            m = rdkit.Chem.rdmolfiles.MolFromMolFile(mol_file)
            if m is None:
                print('\tLoad failed, trying unsanitized structure')
                m = rdkit.Chem.rdmolfiles.MolFromMolFile(mol_file, sanitize=False)
            smiles_str = rdkit.Chem.rdmolfiles.MolToSmiles(m)
            smiles[mol_name] = smiles_str
        
    ''' Optionally output to file '''
    if not smiles_path is None:
        with open(smiles_path, 'w+') as f: 
            f.write('compound\tSMILES\n')
            for mol_file in mol_files:
                mol_name = mol_file.split('/')[-1][:-4]
                if not mol_name in STRUCTURE_BLACKLIST:
                    smiles_str = smiles[mol_name]
                    f.write(mol_name + '\t' + smiles_str + '\n')
    return smiles          

def process_metanetx(mnxm_path=METABOLITES_FILE, 
                     model_path=MODEL_FILE,
                     out_path=METABOLITE_TAG_FILE):
    ''' Extracts the database IDs for each metabolite from MetaNetX '''
    df = pd.read_csv(mnxm_path, header=None, delimiter='\t')
    print('Loaded MetaNetX output:', df.shape)
    model = cobra.io.load_json_model(model_path)
    compartments = model.compartments.keys()
    db_links = {}
    
    for i in range(df.shape[0]):
        if 'MNXM' == df.iloc[i,0][:4]: # metabolite, not biomass
            mnxm, mnxm_name, bigg_all, formula, mw, charge, tags = df.iloc[i,:]
            
            ''' Process the putative BiGG IDs '''
            bigg = bigg_all.split(';')[0][2:] # ignore possible localizations          
            bigg = bigg.replace('_hs_', '_cho_') # replace human label with CHO label
            label = '_'.join(bigg.split('_')[:-1]) # remove "M_" and localization
            db_links[label] = {'mnxm':mnxm, 'bigg':bigg}
                      
            ''' Process IDs to various databases '''
            if type(tags) == str:
                for entry in tags.split(';'):
                    db, tag = entry.split(':')
                    if not db in db_links[label]:
                        db_links[label][db] = tag
                    else:
                        db_links[label][db] += ';' + tag
                        
            ''' Try to pull chemical name from iCHOv1 directly,
                more descriptive than MetaNetX predicted names '''
            name = None
            if bigg in model.metabolites: # if the predicted bigg id matches
                name = model.metabolites.get_by_id(bigg).name
            elif label in NAME_FIXES: # manually curated 
                name = NAME_FIXES[label]
            else: # if need to try all possible bigg id iterations
                candidates = list(map(lambda c: label + '_' + c, compartments))
                for tag in candidates:
                    if tag in model.metabolites:
                        name = model.metabolites.get_by_id(tag).name
                        break;
            if name is None: # still unable to find a good match
                
                print('Unmatched:', bigg_all)
            db_links[label]['name'] = name
            
            ''' Extra information for potential future filtering '''
            db_links[label]['formula'] = formula
            db_links[label]['mw'] = float(mw)
            
    extra_columns = ['name', 'formula', 'mw']
    df_db_links = pd.DataFrame(db_links).transpose()
    df_db_links = df_db_links.reindex(columns=extra_columns + db_priorities)
    df_db_links.to_csv(out_path)
    return df_db_links

def download_mol_files(output_dir=METABOLITE_MOL_DIR, tag_file=METABOLITE_TAG_FILE, 
                      limit=-1, db_order=['chebi','hmdb','kegg','lipidmaps'], 
                      show_all=False, check_dir='../raw_data/CHO_structures_all/'):
    ''' Download structures based on the tag_file. Capable of querying ChEBI
        (TODO: LipidMaps, cannot download directly by manipulating URL) '''
    df = pd.read_csv(tag_file, index_col=0) 
    existing = list(map(lambda x: x[:-4], os.listdir(output_dir)))
    if not check_dir is None:
        existing += list(map(lambda x: x[:-4], os.listdir(check_dir)))
    
    existing = set(existing)
    print('Total Metabolites:', df.shape[0])
    print('    Downloaded:', len(existing))
    print('    Missing:', df.shape[0] - len(existing), '\n')
    counter = 0
    for met in df.index:
        mol_file = output_dir + '/' + met + '.mol'
        mol_file = mol_file.replace('//', '/')
        counter += 1; found = False
        if met in existing:
            if show_all:
                print(met, 'already downloaded')
        else: # otherwise, try to get structure from any of the allowable DBs
            for db in db_order: # test allowable DBs
                if not pd.isnull(df.loc[met,db]): # if a record exists for this DB
                    for tag in df.loc[met,db].split(';'): # test all possible tags
                        ''' Convert the tag into a url'''
                        url = None
                        if db == 'chebi':
                            url = 'https://www.ebi.ac.uk/chebi/saveStructure.do?'
                            url +='defaultImage=true&chebiId=' + str(tag) + '&imageId=0'
                        elif db == 'hmdb':
                            tag = tag.replace('HMDB', 'HMDB00') # URLs pad 0s to ID
                            url = 'http://www.hmdb.ca/structures/metabolites/'
                            url += str(tag) + '.mol'
                        elif db == 'kegg':
                            if tag[0] == 'G': # glycan
                                url = 'https://www.genome.jp/dbget-bin/www_bget?-f+k+glycan+'
                            elif tag[0] == 'C': # generic compound
                                url = 'https://www.genome.jp/dbget-bin/www_bget?-f+m+compound+'
                            url += str(tag)
                        elif db == 'lipidmaps':
                            url = 'https://www.lipidmaps.org/rest/compound/lm_id/'
                            url += str(tag) + '/molfile'
                        
                        ''' Query the url to download the structure '''
                        if not url is None:
                            print(met, '\t', url, end=' ')
                            try:
                                data = urllib.request.urlopen(url, timeout=10)
                                output = ''
                                for bytestring in data:
                                    line = bytestring.decode("utf-8") 
                                    output += line
                                if len(output.strip()) > 0: # got something
                                    if db == 'kegg' and tag[0] == 'G': # KEGG stores glycans as KCF files
                                        out_file = mol_file.replace('.mol','.kcf')
                                    else: 
                                        out_file = mol_file
                                    f = open(out_file, 'w+')
                                    f.write(output)
                                    f.close()
                                    found = True
                                print('')
                            except urllib.error.HTTPError:
                                print('Unable to load URL')

    
                        ''' Stop testing tags if structure was found '''
                        if found: 
                            break
                        
                ''' Stop testing databases if structure was found '''
                if found: 
                    break
                    
        if counter >= limit and limit != -1:
            break;
            
if __name__ == '__main__':
    main()