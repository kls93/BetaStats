
# Parses the AAIndex1 to obtain 566 different continuous, numerical indices
# of amino acid properties, which are then appended to the output DataGen
# dataframe of the CATH / SCOPe domain(s) of interest

import pickle
import pandas as pd
from collections import OrderedDict

aa_positions_dict = OrderedDict({'A': [0, 0],
                                 'R': [1, 0],
                                 'N': [0, 1],
                                 'D': [1, 1],
                                 'C': [0, 2],
                                 'Q': [1, 2],
                                 'E': [0, 3],
                                 'G': [1, 3],
                                 'H': [0, 4],
                                 'I': [1, 4],
                                 'L': [0, 5],
                                 'K': [1, 5],
                                 'M': [0, 6],
                                 'F': [1, 6],
                                 'P': [0, 7],
                                 'S': [1, 7],
                                 'T': [0, 8],
                                 'W': [1, 8],
                                 'Y': [0, 9],
                                 'V': [1, 9]})

# Parses the AAIndex1
with open('docs/AAindex1.txt', 'r') as aa_index_file:
    file_lines = [line for line in aa_index_file.readlines() if line.strip()
                  != '\n']
file_lines = ''.join(file_lines)

aa_indices = OrderedDict()
indices_desc = OrderedDict()
unprocessed_indices = []

index_summaries = file_lines.split('//')
index_summaries = [index_summary for index_summary in index_summaries if
                   index_summary != '\n']
for aa_index in index_summaries:
    id = ''
    desc = ''
    index_dict = ''

    lines = aa_index.split('\n')
    for i, line in enumerate(lines):
        if line[0:1] == 'H':
            id = line.lstrip('H').strip()
        elif line[0:1] == 'D':
            desc = line.lstrip('D').strip()
            count = 1
            while lines[i+count][0:1] != 'R':
                desc += ' ' + lines[i+count].strip()
                count += 1
        elif line[0:1] == 'I':
            index_lines = lines[i+1: i+3]
            index_dict = OrderedDict({'A': '',
                                      'R': '',
                                      'N': '',
                                      'D': '',
                                      'C': '',
                                      'Q': '',
                                      'E': '',
                                      'G': '',
                                      'H': '',
                                      'I': '',
                                      'L': '',
                                      'K': '',
                                      'M': '',
                                      'F': '',
                                      'P': '',
                                      'S': '',
                                      'T': '',
                                      'W': '',
                                      'Y': '',
                                      'V': ''})

            for aa in list(aa_positions_dict.keys()):
                try:
                    index_dict[aa] = float(index_lines[aa_positions_dict[aa][0]].split()[
                                           aa_positions_dict[aa][1]])
                except ValueError:
                    # Deals with NA values
                    index_dict[aa] = index_lines[aa_positions_dict[aa][0]].split()[
                        aa_positions_dict[aa][1]]

    if id != '' and desc != '' and index_dict != '':
        indices_desc[id] = desc
        aa_indices[id] = index_dict
    else:
        unprocessed_indices.append(aa_index)

print('Unprocessed indices: ', unprocessed_indices)

# Loads dataframe output from DataGen and appends AAindex1 values
file_name = input('Specify absolute file path of DataGen dataframe for '
                  'analysis: ')
file_name = ('/Users/ks17361/Lab_work_DW/Beta_structure/Parametrisation/' +
             'Beta_datasets/CATH_2.60.40.10_resn_1.6_rfac_0.2_40_alignment/' +
             'Beta_res_dataframe.pkl')  # DELETE ME!
res_df = pd.read_pickle(file_name)

fasta = res_df['FASTA'].tolist()
for index in list(aa_indices.keys()):
    print('Calculating {} index scores'.format(index))
    index_df = ['']*len(fasta)

    for i, aa in enumerate(fasta):
        if aa in list(aa_positions_dict.keys()):
            index_df[i] = aa_indices[index][aa]
        else:
            index_df[i] = 'NA'

    index_df = pd.DataFrame({'{}'.format(index): index_df})

    res_df = pd.concat([res_df, index_df], axis=1)

res_df.to_csv('Beta_res_dataframe_AA_index.csv')
with open('AAIndices_dicts.pkl', 'wb') as pkl_dicts:
    pickle.dump((aa_indices, indices_desc), pkl_dicts)
