## in this code, we run the scGen for all the tasks
## we can use the below code for running the script
## python scGen_for_all_tasks.py --all_targets task1_ task3_
## the output files will be saved in the ../output/method_outputs/scGen folder
import scanpy as sc
import scgen
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import anndata as ad
from os import listdir
import time
import random
import argparse
import scvi
import os
import re

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

# Create folder to save the results 
save_dirname = '../output/method_outputs/scGen/'
print(save_dirname)
if not os.path.exists(save_dirname): os.makedirs(save_dirname)

def generate_random_seed():
    """
    Generate a random seed.

    Returns:
        int: A random seed.
    """
    return random.randint(0, 2**32 - 1)

def main():
    parser = argparse.ArgumentParser(description='scGen for all tasks')
    parser.add_argument('--all_targets', nargs='+', help='list of target task names')

    args = parser.parse_args()
    all_targets = args.all_targets
    batch_label='batch'
    celltype_label='celltype'

    for target in all_targets:
        print(target)
        ## read data from h5ad files
        pre_dir='../data/process_data/h5ad_files/'      
        onlyfiles = [f for f in listdir(pre_dir) if re.match(target,f)]
        onlyfiles = sorted(onlyfiles)

        # load the data
        adata_a=sc.read_h5ad(pre_dir+onlyfiles[0])
        adata_b=sc.read_h5ad(pre_dir+onlyfiles[1])
        adata_a.obs_names_make_unique()
        adata_b.obs_names_make_unique()
        sc.pp.filter_cells(adata_a, min_genes=200)
        sc.pp.filter_genes(adata_a, min_cells=10)
        sc.pp.filter_cells(adata_b, min_genes=200)
        sc.pp.filter_genes(adata_b, min_cells=10)
        
        # load the one2one orthologous genes mapping table
        if target in ['task4_','task5_','task4-1_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_mouse.csv')
            mapping=mapping.iloc[:,3:5]
            mapping.drop_duplicates(inplace=True,ignore_index=True)   

            # remove all the duplicated genes in mouse and human in the mapping file
            mapping.drop_duplicates(subset=['human_gene'],inplace=True,ignore_index=True)
            mapping.drop_duplicates(subset=['mouse_gene'],inplace=True,ignore_index=True)
            
            ################## convert the mouse gene to human orthologous #############
            human_overlap=list(set(adata_a.var_names) & set(mapping['human_gene']))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['human_gene'])]  
            mapping=mapping.loc[mapping['human_gene'].isin(human_overlap)]  
            adata_b=adata_b[:,adata_b.var_names.isin(mapping['mouse_gene'])]  
            mapping=mapping.loc[mapping['mouse_gene'].isin(adata_b.var_names)]

            mapping=mapping.set_index('mouse_gene')
            mapping=mapping.loc[adata_b.var_names]
            adata_b.var_names=mapping['human_gene']
        elif target in ['task6_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/nema_sty.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            mapping.drop_duplicates(subset=['gene_name'],inplace=True,ignore_index=True)
            mapping.drop_duplicates(subset=['Stylophora'],inplace=True,ignore_index=True)
            nema_overlap=list(set(adata_a.var_names) & set(mapping['gene_name']))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['gene_name'])]  
            mapping=mapping.loc[mapping['gene_name'].isin(nema_overlap)]  
            target_sty_names=adata_b.var_names.str.replace('.{2}$','')
            adata_b.var_names=target_sty_names   

            adata_b=adata_b[:,adata_b.var_names.isin(mapping['Stylophora'])]  
            mapping=mapping.loc[mapping['Stylophora'].isin(adata_b.var_names)]  

            # order the mouse gene order in the mapping file
            mapping=mapping.set_index('Stylophora')
            mapping=mapping.loc[adata_b.var_names]
            adata_b.var_names=mapping['gene_name']

            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]
            
        elif target in ['task8_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/fish_to_frog.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            mapping.drop_duplicates(subset=['fish_gene'],inplace=True,ignore_index=True)
            mapping.drop_duplicates(subset=['frog_gene'],inplace=True,ignore_index=True)

            frog_overlap=list(set(adata_b.var_names) & set(mapping['frog_gene']))
            adata_b=adata_b[:,adata_b.var_names.isin(mapping['frog_gene'])]  
            mapping=mapping.loc[mapping['frog_gene'].isin(frog_overlap)]
            
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['fish_gene'])] 
            
            mapping=mapping.set_index('frog_gene')
            mapping=mapping.loc[adata_b.var_names]
            adata_b.var_names=mapping['fish_gene']

            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        elif target in ['task9_','task9-1_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/ant_to_mouse.csv')
            mapping=mapping.iloc[:,3:5]
            mapping.drop_duplicates(inplace=True,ignore_index=True)

            # remove all the duplicated genes in mouse and human in the mapping file    
            mapping.drop_duplicates(subset=['mouse_gene'],inplace=True,ignore_index=True)
            mapping.drop_duplicates(subset=['ant_gene'],inplace=True,ignore_index=True)
            temp_name=mapping['ant_gene'].str.replace('_','-')
            mapping['ant_gene']=temp_name

            ################## convert the ant gene to mouse orthologous #############
            ant_overlap=list(set(adata_a.var_names) & set(mapping['ant_gene']))
            adata_a=adata_a[:,adata_a.var_names.isin(ant_overlap)]
            mapping=mapping.loc[mapping['ant_gene'].isin(ant_overlap)]  
            mapping=mapping.set_index('ant_gene')
            mapping=mapping.loc[adata_a.var_names]
            
            adata_a.var_names=mapping['mouse_gene']
            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        elif target in ['task10_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/fish_to_fly.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    

            ################## convert the fish gene to fly orthologous #############
            fish_overlap=list(set(adata_a.var_names) & set(mapping['fish_gene']))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['fish_gene'])]  
            mapping=mapping.loc[mapping['fish_gene'].isin(fish_overlap)]  
            fly_overlap=list(set(adata_b.var_names) & set(mapping['fly_symbol']))
            adata_b=adata_b[:,adata_b.var_names.isin(mapping['fly_symbol'])]  
            mapping=mapping.loc[mapping['fly_symbol'].isin(fly_overlap)] 

            mapping=mapping.set_index('fly_symbol')
            mapping=mapping.loc[adata_b.var_names]
            target_orth_fish_names=mapping['fish_gene']
            adata_b.var_names=target_orth_fish_names   
            
            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        elif target in ['task11_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/trichoplax_to_human.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            human_overlap=list(set(adata_a.var_names) & set(mapping['human_gene']))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['human_gene'])] 
            mapping=mapping.loc[mapping['human_gene'].isin(human_overlap)]  
            
            trichoplax_overlap=list(set(adata_b.var_names) & set(mapping['trichoplax_gene']))
            adata_b=adata_b[:,adata_b.var_names.isin(mapping['trichoplax_gene'])] 
            mapping=mapping.loc[mapping['trichoplax_gene'].isin(trichoplax_overlap)]  
            mapping=mapping.set_index('trichoplax_gene')
            mapping=mapping.loc[adata_b.var_names]
            target_orth_names=mapping['human_gene']
            adata_b.var_names=target_orth_names   
            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        elif target in ['task12_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/ciona_to_nema.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            ciona_overlap=list(set(adata_a.var_names) & set(mapping['ciona_gene']))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['ciona_gene'])] 
            mapping=mapping.loc[mapping['ciona_gene'].isin(ciona_overlap)]   
            
            nema_overlap=list(set(adata_b.var_names) & set(mapping['nema_gene']))
            adata_b=adata_b[:,adata_b.var_names.isin(mapping['nema_gene'])] 
            mapping=mapping.loc[mapping['nema_gene'].isin(nema_overlap)] 
            mapping=mapping.set_index('nema_gene')
            mapping=mapping.loc[adata_b.var_names]
            target_orth_names=mapping['ciona_gene']

            adata_b.var_names=target_orth_names   
            adata_b.var_names_make_unique()
            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        elif target in ['task13_','task33_','task34_','task35_','task36_','task37_','task38_','task39_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/sea_urchin_to_zebrafish.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            mapping = mapping.drop_duplicates(subset="sea_urchin_gene")
            mapping = mapping.drop_duplicates(subset="fish_gene")
            a_overlap=list(set(adata_a.var_names) & set(mapping['sea_urchin_gene']))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['sea_urchin_gene'])]  
            mapping=mapping.loc[mapping['sea_urchin_gene'].isin(a_overlap)]  

            b_overlap=list(set(adata_b.var_names) & set(mapping['fish_gene']))
            adata_b=adata_b[:,adata_b.var_names.isin(mapping['fish_gene'])] 
            mapping=mapping.loc[mapping['fish_gene'].isin(b_overlap)] 
            mapping=mapping.set_index('fish_gene')
            mapping=mapping.loc[adata_b.var_names]
            target_orth_names=mapping['sea_urchin_gene']

            adata_b.var_names=target_orth_names   
            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        elif target in ['task14_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/celegan_to_human.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            mapping = mapping.drop_duplicates(subset="human_gene")
            a_overlap=list(set(adata_a.var_names) & set(mapping['celegan_gene']))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['celegan_gene'])]  
            mapping=mapping.loc[mapping['celegan_gene'].isin(a_overlap)] 
            
            b_overlap=list(set(adata_b.var_names) & set(mapping['human_gene']))
            adata_b=adata_b[:,adata_b.var_names.isin(mapping['human_gene'])]  
            mapping=mapping.loc[mapping['human_gene'].isin(b_overlap)] 
            mapping=mapping.set_index('human_gene')
            mapping=mapping.loc[adata_b.var_names]
            target_orth_names=mapping['celegan_gene']

            adata_b.var_names=target_orth_names   
            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        elif target in ['task15_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/octopus_to_human.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            mapping = mapping.drop_duplicates(subset="human_gene")
            a_overlap=list(set(adata_a.var_names) & set(mapping['human_gene']))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['human_gene'])] 
            mapping=mapping.loc[mapping['human_gene'].isin(a_overlap)]  
            
            b_overlap=list(set(adata_b.var_names) & set(mapping['octopus_gene']))
            adata_b=adata_b[:,adata_b.var_names.isin(mapping['octopus_gene'])]  
            mapping=mapping.loc[mapping['octopus_gene'].isin(b_overlap)]  
            mapping=mapping.set_index('octopus_gene')
            mapping=mapping.loc[adata_b.var_names]
            target_orth_names=mapping['human_gene']

            adata_b.var_names=target_orth_names   
            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        elif target in ['task16_']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/schmid_to_human.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            mapping = mapping.drop_duplicates(subset="human_gene")
            a_overlap=list(set(adata_a.var_names) & set(mapping['human_gene']))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping['human_gene'])] 
            mapping=mapping.loc[mapping['human_gene'].isin(a_overlap)]  
            
            b_overlap=list(set(adata_b.var_names) & set(mapping['schmidtea_gene']))
            adata_b=adata_b[:,adata_b.var_names.isin(mapping['schmidtea_gene'])]  
            mapping=mapping.loc[mapping['schmidtea_gene'].isin(b_overlap)] 
            mapping=mapping.set_index('schmidtea_gene')
            mapping=mapping.loc[adata_b.var_names]
            target_orth_names=mapping['human_gene']

            adata_b.var_names=target_orth_names   
            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        elif target in ["task" + str(i) + "_" for i in range(17,33)]:
            if target in ['task17_']:
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/cat_to_tiger.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)    
                a_name='cat_gene'
                b_name='tiger_gene'
            elif target in ['task18_']:
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/cat_to_dog.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)    
                a_name='cat_gene'
                b_name='dog_gene'
            elif target == 'task19_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_to_MF.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)    
                a_name='MF_gene'
                b_name='human_gene'
            elif target == 'task20_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_to_MM.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)    
                a_name='MM_gene'
                b_name='human_gene'
            elif target == 'task21_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_to_mouse_mm10.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)    
                a_name='human_gene'
                b_name='mouse_gene'
            elif target == 'task22_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/mouse_to_MF.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='MF_gene'
                b_name='mouse_gene'
            elif target == 'task23_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/mouse_to_MM.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='MM_gene'
                b_name='mouse_gene'
            elif target == 'task24_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_to_pig.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='human_gene'
                b_name='pig_gene'
            elif target == 'task25_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/mouse_to_pig.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='mouse_gene'
                b_name='pig_gene'
            elif target == 'task26_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/pig_to_MF.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='MF_gene'
                b_name='pig_gene'
            elif target == 'task27_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/pig_to_MM.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='MM_gene'
                b_name='pig_gene'
            elif target == 'task28_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_to_fish.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='fish_gene'
                b_name='human_gene'
            elif target == 'task29_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/fish_to_MF.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='MF_gene'
                b_name='fish_gene'
            elif target == 'task30_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/fish_to_cat.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='cat_gene'
                b_name='fish_gene'
            elif target == 'task31_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/pig_to_fish.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='fish_gene'
                b_name='pig_gene'
            elif target == 'task32_':
                mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/fish_to_MM.csv')
                mapping.drop_duplicates(inplace=True,ignore_index=True)
                a_name='MM_gene'
                b_name='fish_gene'

            mapping = mapping.drop_duplicates(subset=a_name)
            a_overlap=list(set(adata_a.var_names) & set(mapping[a_name]))
            adata_a=adata_a[:,adata_a.var_names.isin(mapping[a_name])] 
            mapping=mapping.loc[mapping[a_name].isin(a_overlap)]  
            
            b_overlap=list(set(adata_b.var_names) & set(mapping[b_name]))
            adata_b=adata_b[:,adata_b.var_names.isin(mapping[b_name])]  
            mapping=mapping.loc[mapping[b_name].isin(b_overlap)]
            mapping = mapping.drop_duplicates(subset=b_name) 
            mapping=mapping.set_index(b_name)
            mapping=mapping.loc[adata_b.var_names]
            target_orth_names=mapping[a_name]
            adata_b.var_names=target_orth_names   
            common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
            adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
            adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]

        adata_a.layers['counts'] = adata_a.X.copy()
        adata_b.layers['counts'] = adata_b.X.copy()


        adata_ls = ad.concat([adata_a,adata_b],join='inner')
        
        adata_ls.obs_names_make_unique()

        random.seed(time.time())
        seed = generate_random_seed()

        random.seed(seed)
        np.random.seed(seed)
        scvi.settings.seed = seed
        start_time=time.time()

        sc.pp.normalize_total(adata_ls)
        sc.pp.log1p(adata_ls)
        sc.tl.pca(adata_ls,svd_solver='arpack')
        
        adata_ls.obs['cell_type']=adata_ls.obs[celltype_label].tolist()
        adata_ls.obs['batch']=adata_ls.obs[batch_label].tolist()
        adata_ls.obs['batch'] = adata_ls.obs['batch'].astype(str)

        scgen.SCGEN.setup_anndata(adata_ls, batch_key="batch", labels_key="cell_type")
        
        # create model and save it
        model = scgen.SCGEN(adata_ls)
        model.save(os.path.join(save_dirname,re.sub('_','',target)+"_seed_"+str(seed)+'.pt'), overwrite=True)
    
        # train model
        model.train(
            max_epochs=100,
            batch_size=32,
            early_stopping=True,
            early_stopping_patience=25)
        
        # batch removal
        corrected_adata = model.batch_removal()
        #corrected_adata
        end_time=time.time()
        elapsed_time = end_time - start_time
        with open(f"../output/method_outputs/scGen/scGen_timing_log.txt", "a") as file:
            file.write(f"{re.sub('_', '', target)}, seed: {seed}, Elapsed time: {elapsed_time:.2f} seconds\n")

        sc.pp.neighbors(corrected_adata, use_rep="corrected_latent")
        corrected_adata.obs['orig.ident']=corrected_adata.obs['batch']
        corrected_adata.write_h5ad(os.path.join(save_dirname,target+str(seed)+'.h5ad'))


        # save the pca files
        filename = target+str(seed)+'_scGen_embedding.txt'
        coln_pca = []
        for i in range(corrected_adata.obsm['corrected_latent'].shape[1]):
            coln_pca.append("X_scgen"+str(i+1))
        coln_pca.append('batchlb')
        coln_pca.append('celltype')
        target_obj=pd.DataFrame(corrected_adata.obsm['corrected_latent'])
        target_obj['batch']=corrected_adata.obs['batch'].values
        target_obj['celltype']=corrected_adata.obs['cell_type'].values
        target_obj.columns=coln_pca
        adata_ls.obs_names_make_unique()
        target_obj.index=adata_ls.obs_names
        target_obj.to_csv(save_dirname+filename)

main()
