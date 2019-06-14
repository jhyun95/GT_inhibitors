from chemspipy import ChemSpider
cs = ChemSpider('<YOUR-API-KEY>')
import pandas as pd

def get_SMILES(compounds):
    """ Query ChemSpider for SMILES """
    smiles = []
    for compound in compounds:
        try:
            float(compound)
            smiles.append(compound)
        except:
            result = cs.search(compound)
            smiles.append([i.smiles for i in result])
    return smiles

if __name__ == '__main__':
    """ Load GT inhibitors/activators and query ChemSpider for SMILES """
    inhib_df = pd.read_csv('../data/inhibitors.csv')
    activ_df = pd.read_csv('../data/activators.csv')
    inhib_smiles_df = pd.DataFrame(columns=inhib_df.columns)
    for column in inhib_df.columns:
        inhib_smiles_df[column] = get_SMILES(inhib_df[column])
    activ_smiles_df = pd.DataFrame(columns=activ_df.columns)
    for column in activ_df.columns:
        activ_smiles_df[column] = get_SMILES(activ_df[column])
    inhib_smiles_df.to_csv('../data/inhibitors_smiles.csv', index=False)
    activ_smiles_df.to_csv('../data/activators_smiles.csv', index=False)