
import os
import shutil
import pickle
import PIL
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# draws plots of every amino acid index
# vs. 3 structural features (buried surface area (and core vs. surface),
# residue displacement from the centre of the strand (both absolute and %, as
# both a discrete and a continuous variable), and edge vs. central)
res_df = pd.read_csv('Beta_res_dataframe_AA_index.csv')
with open('AAIndices_dicts.pkl', 'rb') as pkl_dicts:
    aa_indices, indices_desc = pickle.load(pkl_dicts)

if os.path.isdir('Figures'):
    shutil.rmtree('Figures')
os.mkdir('Figures')

structural_features = ['EDGE_OR_CNTRL', 'BURIED_SURFACE_AREA()%)',
                       'STRAND_POS(ABS)', 'STRAND_POS(%)']
colours = ['b', 'g', 'r', 'y']
count = -1

for feature in structural_features:
    count += 1
    for index in list(aa_indices.keys()):
        print('Generating scatter plot of {} vs. {}'.format(feature, index))

        x = res_df['{}'.format(index)].tolist()
        y = res_df['{}'.format(feature)].tolist()

        if feature in ['EDGE_OR_CNTRL']:
            new_y = []
            for entry in y:
                if entry in ['edge', 'core']:
                    new_y.append(1)
                elif entry in ['central', 'surface']:
                    new_y.append(0)
                else:
                    new_y.append('NA')
            y = new_y

        plot_df = pd.DataFrame({'{}'.format(index): x,
                                '{}'.format(feature): y})
        plot_df = plot_df[~plot_df['{}'.format(index)].isin(['', 'NA'])]
        plot_df = plot_df.reset_index(drop=True)
        plot_df = plot_df[~plot_df['{}'.format(feature)].isin(['', 'NA'])]
        plot_df = plot_df.reset_index(drop=True)

        plt.clf()
        sns.regplot(x=plot_df['{}'.format(index)],
                    y=plot_df['{}'.format(feature)], ci=None, marker='.',
                    scatter_kws={'alpha': 0.05}, color=colours[count])
        plt.xlabel('{}'.format(indices_desc[index]))
        plt.ylabel('{}'.format(feature.lower()))
        plt.margins(0.02)
        plt.savefig('Figures/{}_vs_{}.png'.format(feature.lower(), index))
