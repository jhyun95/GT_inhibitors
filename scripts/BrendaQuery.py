import hashlib
from SOAPpy import SOAPProxy
import re
import pandas as pd

def extract_ECs(EC_list_file):
    relevant_ECs = {}
    with open(EC_list_file, 'r+') as f:
        for line in f:
            enzyme, EC = line.strip().split(' - ')
            relevant_ECs[enzyme] = EC
    return relevant_ECs
def search_Inhibitors_and_Activators(EC):
    endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
    password = hashlib.sha256("PASSWORD").hexdigest()
    parameters = "EMAIL,"+password+",ecNumber*"+EC
    client = SOAPProxy(endpointURL)
    resultString = client.getInhibitors(parameters)
    inhibitors = set([match.group(1) for match in re.finditer('#inhibitor\*(.+?)#', resultString)])
    resultString = client.getActivatingCompound(parameters)
    activators = set([match.group(1) for match in re.finditer('#activatingCompound\*(.+?)#', resultString)])
    return inhibitors, activators

if __name__ == '__main__':
    ECs = extract_ECs('../data/GT_ECs.txt')
    d = {}
    g = {}
    for ec in ECs.values():
        inhibitors, activators = search_Inhibitors_and_Activators(ec)
        d[ec] = list(inhibitors)
        g[ec] = list(activators)
    inhib_df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in d.iteritems()]))
    activ_df = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in g.iteritems()]))
    inhib_df.to_csv('../data/inhibitors.csv', index=False)
    activ_df.to_csv('../data/activators.csv', index=False)
