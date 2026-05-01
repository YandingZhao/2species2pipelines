## this file to gather all the evaluation results together for each method
## the output file will be saved in the ../output/final_evaluation_table folder

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
from os import listdir
from os.path import isfile, join
import re
import os

save_dir = '../output/final_evaluation_table/'
print(save_dir)
if not os.path.exists(save_dir): os.makedirs(save_dir)

pre_dir='../output/evaluation/'

all_targets=["task" + str(i) for i in range(3,39)]
all_targets.insert(7,'task9-1')

unique_celltype_number={'task4': 72,'task3': 12,'task11': 12,'task9': 2,'task16': 17,
'task10': 5,'task8': 50,'task8-2': 50,'task9-1': 2,'task13': 12,'task7': 11,'task15': 12,'task6': 13,'task12': 13,
'task14': 16,'task17': 3, 'task17-2': 3,'task18': 3, 'task18-2': 3,'task19': 10, 'task20':7, 'task21':11, 'task22':11,'task23':8,'task24':9,'task25':7,
'task26':8,'task27':7,'task28':2,'task29':2,'task30':3,'task31':2,'task32':2,'task33':12,'task34':12,'task35':12,'task36':12,'task37':12,'task38':12,'task39':12}

unique_batch_number={'task4': 11,'task3': 2,'task11': 2,'task9': 2,'task16': 2,'task10': 2,
'task8': 2,'task8-2': 2,'task9-1': 2,'task13': 2,'task7': 2,'task15': 2,'task6': 2,'task12': 2,'task14': 2, 'task17': 2, 'task17-2': 2,
'task18': 2, 'task18-2': 2,'task19': 2, 'task20':2, 'task21':2, 'task22':2,'task23':2,'task24':2,'task25':2,'task26':2,'task27':2,'task28':2,
'task29':2,'task30':2,'task31':2,'task32':2,'task33':2,'task34':2,'task35':2,'task36':2,'task37':2,'task38':2,'task39':2}

#all_methods=listdir(pre_dir)
all_methods=['bbknn','fastmnn','harmony','scVI','seurat4','scanorama','scGen','samap']
#all_methods=['saturn']
for method in all_methods:
    print(method)
    final_target={'ari_batch':[],'ari_celltype':[],'nmi_batch':[],'nmi_celltype':[],
                'asw_batch':[],'asw_celltype':[],'lisi_batch':[],'lisi_celltype':[],
                'kBET':[],'PCR_batch':[],'isolated_label_F1':[],'isolated_label_silhouette':[],'graph_conn':[],'hvg_overlap':[],'trajectory':[]}


    for target in all_targets:
        print(target)
        # get the ARI result
        ari=pd.read_table(pre_dir+method+'/'+target+'_'+method+'_ARI.txt')
        # get the last row in ari which is the mean value
        ari_mean=ari.iloc[-1]
        if ari_mean['ari_batch'] < 0:
            ari_mean['ari_batch']=0
 
        # append the adjusted rank value to the final_target 
        final_target['ari_batch'].append((1-ari_mean['ari_batch']))
        final_target['ari_celltype'].append(ari_mean['ari_celltype'])

        # get the NMI result which are already normalized
        nmi=pd.read_csv(pre_dir+method+'/'+target+'_'+method+'_NMI.csv')
        final_target['nmi_batch'].append(1-nmi['nmi_value_batch'].iloc[0])
        final_target['nmi_celltype'].append(nmi['nmi_value_celltype'].iloc[0])
    
        # get the ASW result which are already normalized
        if method not in ['samap','bbknn']:
            asw=pd.read_csv(pre_dir+method+'/'+target+'_'+method+'_ASW_metric.csv')
            final_target['asw_batch'].append(asw['asw_batch_norm_sub'].iloc[-1])
            final_target['asw_celltype'].append(asw['asw_celltype_norm'].iloc[-1])
        else:
            final_target['asw_batch'].append(np.nan)
            final_target['asw_celltype'].append(np.nan)
        
        # get the lisi result
        if method not in ['bbknn','samap']:
            lisi_batch_df=pd.read_table(pre_dir+method+'/'+target+'_'+method+'_lisi_batch_40.txt')
            lisi_celltype_df=pd.read_table(pre_dir+method+'/'+target+'_'+method+'_lisi_celltype_40.txt')

            # normalize the dataframe by the unique number of cell type and batch
            lisi_batch_df=lisi_batch_df.iloc[:-1]
            lisi_celltype_df=lisi_celltype_df.iloc[:-1]
            lisi_batch_df['batch']=(lisi_batch_df['batch']-1)/(unique_batch_number[target]-1)
            lisi_celltype_df['cell_type']=(unique_celltype_number[target]-lisi_celltype_df['cell_type'])/(unique_celltype_number[target]-1)
            # get the mean lisi result:
            lisi_batch_mean=lisi_batch_df['batch'].mean()
            lisi_celltype_mean=lisi_celltype_df['cell_type'].mean()

            final_target['lisi_batch'].append(lisi_batch_mean)
            final_target['lisi_celltype'].append(lisi_celltype_mean)  #should not be substracted by 1, the normalization form already did it
        else:
            lisi=pd.read_csv(pre_dir+method+'/'+target+'_'+method+'_LISI.csv')
            lisi=lisi.iloc[:-1]
            batch_column_index = lisi.columns[1]
            celltype_column_index = lisi.columns[2]
            lisi_batch_df=lisi[[batch_column_index]]
            lisi_celltype_df=lisi[[celltype_column_index]]
            lisi_batch_df['batch']=(lisi_batch_df[batch_column_index]-1)/(unique_batch_number[target]-1)
            lisi_celltype_df['cell_type']=(unique_celltype_number[target]-lisi_celltype_df[celltype_column_index])/(unique_celltype_number[target]-1)
            lisi_batch_mean=lisi_batch_df['batch'].mean()
            lisi_celltype_mean=lisi_celltype_df['cell_type'].mean()
            final_target['lisi_batch'].append(lisi_batch_mean)
            final_target['lisi_celltype'].append(lisi_celltype_mean)
            
        # get the result from scib:kBet, PCR_batch, isolated_F1, isolated_silhouette, graph_connectivity, hvg_overlap, trajectory, true_tranjrctory
        scib_output=pd.read_csv(pre_dir+method+'/'+target+'_'+method+'_scib_output.csv',index_col=0)
        final_target['kBET'].append(scib_output.loc['kBET'][0])
        final_target['PCR_batch'].append(scib_output.loc['PCR_batch'][0])
        final_target['graph_conn'].append(scib_output.loc['graph_conn'][0])
        final_target['hvg_overlap'].append(scib_output.loc['hvg_overlap'][0])
        final_target['trajectory'].append(scib_output.loc['trajectory'][0])
    evaluation_df=pd.DataFrame.from_dict(final_target)
    evaluation_df.index=all_targets
    evaluation_df.to_csv(save_dir+method+'_final_task3_38.csv')
