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
        f_out = open(out_EC_data, 'w+', encoding='utf8')
        recording = False
        for line in f_in:
            if recording: # in the middle of relevant EC
                f_out.write(line)
                recording = not(line.strip() == '///')
            else: # outside a relevant EC
                if line[:3] == 'ID\t': # check if entering a new EC
                    if line.split()[1] in relevant_ECs:
                        relevant_ECs.remove(line.split()[1])
                        f_out.write(line)
                        recording = True
                        print(line.strip())
        f_out.close()
    print('Unidentified ECs:', relevant_ECs)
    
def extract_inhibitors(EC_data='../data/brenda_GTs.txt',
                       out_table='../data/inhibitors.tsv'):
    ''' Extracts inhibitors from raw BRENDA data, along with their
        target EC and KI data/comments if available '''
        
    def process_inhibitor_text(data):
        ''' Extracts the name and commentary from raw BRENDA text '''
        full_line = '\n'.join(data)
        full_line = full_line[3:] # drop the flag and tab separator
        full_line = full_line[full_line.index('#')+1:] # drop protein ref tag
        full_line = full_line[full_line.index('#')+1:] # drop protein ref tag
        if '(#' in full_line and '>)' in full_line: # comment is present
            comment_start = full_line.index('(#')
            comment_end = full_line.index('>)')
            comment = full_line[comment_start+1:comment_end+1]
            inhibitor_name = full_line[:comment_start].strip()
        else:
            comment = ''
            inhibitor_name = full_line[:full_line.index('<')].strip()
        #TODO: line breaks appear arbitrary? commas between numbers, spaces between text/number
        comment = comment.replace('\n',' ')
        inhibitor_name = inhibitor_name.replace('\n',',')
        return (inhibitor_name, comment)
                                              
    entries = [] # each entry is (inhibitor name, EC, KI, comments)
    current_EC = None; current_section = None; current_data = []
    with open(EC_data, 'r+', encoding='utf8') as f:
        for line in f:
            if current_section == 'IN': # in the middle of an inhibitor entry
                if len(line.strip()) == 0 or line[:2] == 'IN': # entry ended
                    print('----------------------------------------')
                    print(current_data)
                    print(process_inhibitor_text(current_data))
                    current_section = None    
                else: # continuing entry
                    current_data.append(line.strip())
                
            elif current_section == 'KI': # in the middle of a KI entry
                if len(line.strip()) == 0 or line[:2] == 'KI': # entry ended
                    # TODO: parsing KI values
                    current_section = None
                else: # continuing entry
                    current_data.append(line.strip()) 
            
            if current_section is None and len(line) > 2: # between sections
                key = line[:3] # line label
                if key == 'ID\t': # entering a new EC section
                    current_EC = line.split()[1]
                elif key == 'IN\t': # entering new inhibitor section
                    current_section = 'IN'
                    current_data = [line.strip()]
                elif key == 'KI\t': # entering new KI data section
                    current_section = 'KI'
                    current_data = [line.strip()]
        
if __name__ == '__main__':  
#    extract_ECs()
    extract_inhibitors()