
# Parses the AAIndex(1) to obtain 566 different continuous, numerical indices of
# amino acid properties, which are then plotted against the structural features
# of interest

import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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

# Parses the AAIndex(1)
with open('docs/AAindex1.txt', 'r') as aa_index_file:
    file_lines = [line for line in aa_index_file.readlines() if line.strip()
                  != '\n']
file_lines = ''.join(file_lines)

aa_indices = OrderedDict()
indices_desc = OrderedDict()
unprocessed_indices = []

index_summaries = file_lines.split('//')
for aa_index in index_summaries:
    id = ''
    desc = ''
    index_dict = ''

    lines = aa_index.split('\n')
    for i, line in enumerate(lines):
        if line.startswith('H'):
            id = line.strip('H').strip()
        elif line.startswith('D'):
            desc = line.strip('D').strip()
        elif line.startswith('I'):
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

# Loads dataframe output from DataGen and draws plots of every amino acid index
# vs. 3 structural features (buried surface area (and core vs. surface),
# residue displacement from the centre of the strand (both absolute and %, as
# both a discrete and a continuous variable), and edge vs. central)
if os.path.isdir('Figures'):
    shutil.rmtree('Figures')
os.mkdir('Figures')

res_df = pd.read_pickle(
    '/Users/ks17361/Lab_work_DW/Beta_structure/Parametrisation/Beta_datasets/CATH_2.60.40.10_resn_1.6_rfac_0.2_40_alignment/Beta_res_dataframe.pkl'
)
fasta = res_df['FASTA'].tolist()
for index in list(aa_indices.keys()):
    print('Calculating {} index scores'.format(index))
    index_df = ['']*len(fasta)

    for i, aa in enumerate(fasta):
        if aa in list(aa_positions_dict.keys()):
            index_df[i] = aa_indices[index][aa]

    index_df = pd.DataFrame({'{}'.format(index): index_df})

    res_df = pd.concat([res_df, index_df], axis=1)

res_df.to_csv('Beta_res_dataframe_AA_index.csv')

structural_features = ['EDGE_OR_CNTRL', 'BURIED_SURFACE_AREA()%)',
                       'STRAND_POS(ABS)', 'STRAND_POS(%)', 'CORE_OR_SURFACE']
for index in list(aa_indices.keys()):
    for feature in structural_features:
        x = res_df['{}'.format(index)]
        y = res_df['{}'.format(feature)]

        if feature in ['EDGE_OR_CNTRL', 'CORE_OR_SURFACE']:
            new_y = []
            for entry in y:
                if entry in ['edge', 'core']:
                    new_y.append(1)
                elif entry in ['central', 'surface']:
                    new_y.append(0)
                else:
                    new_y.append('NA')
        y = new_y
        df = pd.DataFrame({'{}'.format(index): x,
                           '{}'.format(feature): y})
        df = df[~df['{}'.format(index)].isin(['', 'NA'])]
        df = df.reset_index(drop=True)
        df = df[~df['{}'.format(feature)].isin(['', 'NA'])]
        df = df.reset_index(drop=True)

        plt.clf()
        plt.plot(x=df['{}'.format(index)].tolist(), y=df['{}'.format(feature)].tolist())
        plt.show()
        import sys
        sys.exit()
        plt.savefig()
