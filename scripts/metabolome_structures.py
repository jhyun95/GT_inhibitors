#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 14:47:27 2019

@author: jhyun95
"""

import os, urllib
import cobra
import pandas as pd

MODEL_FILE = '../data/iCHOv1.json'
METABOLITES_FILE = '../data/iCHOv1.chemicals.tsv'
METABOLITE_TAG_FILE = '../data/iCHOv1.chemical_tags.csv'
OUTPUT_DIR = '../raw_data/CHO_structures/'

NAME_FIXES = { # for BiGG IDs that are on MetaNetX but not iCHOv1?
    'dhnpt': '7,8-dihydroneopterin',
    'sdhlam': 'S-succinyldihydrolipoamide',
    'ficytc': 'ferricytochrome',
    'pepd': 'peptide',
    'm2macchitppdol': '(alpha-D-mannosyl)2-beta-D-mannosyldiacetylchitobiosyldiphosphodolichol',
    'm1macchitppdol': 'alpha-D-mannosyl-beta-D-mannosyl-diacylchitobiosyldiphosphodolichol',
    'm3macchitppdol': '(alpha-D-mannosyl)3-beta-D-mannosyl-diacetylchitodiphosphodolichol',
    'ptd2meeta_SC': 'phosphatidyldi-N-methylethanolamine_tomerge'}

db_priorities = ['chebi', 'hmdb', 'kegg', 'lipidmaps'] # structures readily available
db_priorities += ['seed', 'metacyc', 'sabiork', 'reactome', 'mnxm', 'bigg'] # structures less available

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

def download_mol_files(output_dir=OUTPUT_DIR, tag_file=METABOLITE_TAG_FILE, 
                      limit=-1, db_order=['chebi'], show_all=False):
    ''' Download structures based on the tag_file. Capable of querying ChEBI
        (TODO: HMDB, KEGG, and LipidMaps querying) '''
    df = pd.read_csv(tag_file, index_col=0) 
    existing = list(map(lambda x: x[:-4], os.listdir(output_dir)))
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
                            pass # TODO
                        elif db == 'kegg':
                            pass # TODO
                        elif db == 'lipidmaps':
                            pass # TODO
                        
                        ''' Query the url to download the structure '''
                        if not url is None:
                            print(met, '\t', url)
                            
                            data = urllib.request.urlopen(url, timeout=10)
                            output = ''
                            for bytestring in data:
                                line = bytestring.decode("utf-8") 
                                output += line
                            if len(output.strip()) > 0: # got something
                                f = open(mol_file, 'w+')
                                f.write(output)
                                f.close()
                                found = True
    
                        ''' Stop testing tags if structure was found '''
                        if found: 
                            break
                        
                ''' Stop testing databases if structure was found '''
                if found: 
                    break
                    
        if counter >= limit and limit != -1:
            break;
            

#df = process_metanetx()
download_mol_files()
