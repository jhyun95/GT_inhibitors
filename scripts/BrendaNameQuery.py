import hashlib
from SOAPpy import SOAPProxy
import re
import pandas as pd
import pickle

def extract_ECs(EC_list_file):
    relevant_ECs = {}
    with open(EC_list_file, 'r+') as f:
        for line in f:
            enzyme, EC = line.strip().split(' - ')
            relevant_ECs[enzyme] = EC
    return relevant_ECs
def search_Inhibitor_names(EC):
    endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
    password = hashlib.sha256("PASSWORD").hexdigest()
    parameters = "EMAIL,"+password+",ecNumber*"+EC
    client = SOAPProxy(endpointURL)
    resultString = client.getInhibitors(parameters)
    inhibitors = [match.group(1) for match in re.finditer('#inhibitor\*(.*?)#', resultString)]
    inhibitor_structureIDs = [match.group(1) for match in re.finditer('#ligandStructureId\*(\d+)#', resultString)]
    d = dict(zip(inhibitor_structureIDs, inhibitors))
    return d

if __name__ == '__main__':
    inhib_name_df = pd.DataFrame(columns=['EC', 'group_id', 'name'])
    ECs = extract_ECs('../data/GT_ECs.txt')
    i = 0
    for ec in ECs.values():
        d = search_Inhibitor_names(ec)
        for k,v in d.iteritems():
            inhib_name_df.loc[i] = [ec, k, v]
            i += 1
    inhib_name_df.to_csv('../data/inhibitor_names.csv', index=False)