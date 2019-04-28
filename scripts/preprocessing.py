# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 18:04:07 2019

@author: jhyun_000
"""

def extract_ECs(EC_list_file='../data/GT_ECs.txt', 
                full_EC_data='../raw_data/brenda_download.txt',
                out_EC_data='../data/brenda_GTs.txt'):
    ''' Reduces full BRENDA data dump to that of the relevant ECs '''
    relevant_ECs = []
    with open(EC_list_file, 'r+') as f:
        for line in f:
            enzyme, EC = line.strip().split(' - ')
            relevant_ECs.append(EC)
    print('Loaded ECs:', relevant_ECs)
    
    with open(full_EC_data, 'r+', encoding='utf8') as f_in:
        f_out = open(out_EC_data, 'w+')
        recording = False
        for line in f_in:
            if recording: # in the middle of relevant EC
                f_out.write(line)
                recording = not(line.strip() == '///')
            else: # outside a relevant EC
                if line[:2] == 'ID': # check if entering a new EC
                    if line.split()[1] in relevant_ECs:
                        relevant_ECs.remove(line.split()[1])
                        f_out.write(line)
                        recording = True
                        print(line.strip())
        f_out.close()
    print('Unidentified ECs:', relevant_ECs)
    
extract_ECs()