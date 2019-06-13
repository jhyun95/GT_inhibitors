import pandas as pd
import subprocess
import networkx as nx

if __name__ == '__main__':
    inhibitor_smiles = pd.read_csv('../data/inhibitor_smiles.tsv', sep='\t')
    inhibitor_smiles_processed = inhibitor_smiles[['SMILES', 'group_id']]
    inhibitor_smiles_processed.columns = ['SMILES', 'compound']
    cho_smiles = pd.read_csv('../data/iCHOv1_smiles.tsv', sep='\t')
    cho_smiles_processed = cho_smiles[['SMILES', 'compound']]
    df = pd.concat([inhibitor_smiles_processed, cho_smiles_processed])
    
    df.to_csv('../data/temp.smi', sep='\t', header=False, index=False)
    df_generate = subprocess.run(['sng', 'generate', '-p', '100', '-o', '../data/temp.tmp', '../data/temp.smi'], check=True, text=True)
    df_aggregate = subprocess.run(['sng', 'aggregate', '-m', '../data/temp_m.tmp', '-o', '../data/temp.network', '../data/temp.tmp'], check=True, text=True)
    network_df = pd.read_csv('../data/temp.network', sep='\t')
    molecule_map = pd.read_csv('../data/temp_m.tmp', sep='\t')
    G = nx.DiGraph()
    G.add_nodes_from([str(i) for i in list(network_df['ID'])])
    for index in network_df.index:
        source = str(network_df['ID'][index])
        try:
            targets = [i for i in network_df['SUBSCAFFOLDS'][index].split(',') if len(i)]
            for target in targets:
                G.add_edges_from([(source, str(target))])
        except: 
            pass
    molecules = set(molecule_map['MOLECULE_ID'])
    scaffold_d = dict(zip(set(network_df['SMILES']), [0]*len(set(network_df['SMILES']))))
    molecule_scaffold_d = dict(zip(molecules,[scaffold_d.copy() for i in range(len(molecules))]))
    for index in molecule_map.index:
        molecule = molecule_map['MOLECULE_ID'][index]
        source = str(molecule_map['SCAFFOLD_ID'][index])
        descendants =  nx.descendants(G, source)
        molecule_scaffold_d[molecule][list(network_df.loc[network_df['ID'] == int(source)]['SMILES'])[0]] = 1
        for descendant in descendants:
            molecule_scaffold_d[molecule][list(network_df.loc[network_df['ID'] == int(descendant)]['SMILES'])[0]] = 1
    removal = subprocess.run(['rm', '-r', '../data/temp.smi', '../data/temp.tmp', '../data/temp.network', '../data/temp_m.tmp'], check=True, text=True)
    molecule_scaffold_df = pd.DataFrame(data=molecule_scaffold_d).T
    temp = molecule_scaffold_df[molecule_scaffold_df.index.isin([str(i) for i in inhibitor_smiles_processed['compound']])].T
    temp.to_csv('../data/scaffold_fingerprints_inhibitors_CHO_i.csv')
    temp = molecule_scaffold_df[molecule_scaffold_df.index.isin([str(i) for i in cho_smiles_processed['compound']])].T
    temp.to_csv('../data/scaffold_fingerprints_inhibitors_CHO_m.csv')