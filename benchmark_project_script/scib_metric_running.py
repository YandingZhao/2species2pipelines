# use the benchmark_py environment
# below code is to run the evaluation on kBET, PCR batch, graph connectivity, Highly variable gene conservation, trajectory conservation using scib package
# for usage of the code:
# python scib_metric_running.py --method <method_name> --target <task_id> 
import os
import sys
import copy
import glob
import anndata as ad
import scib
import argparse
import pandas as pd
import numpy as np
import re
import scanpy as sc

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='scib evaluation')
    parser.add_argument('--method',type=str,help='the method name')
    parser.add_argument('--target',type=str,help='the task id')
    args = parser.parse_args()

    method=args.method
    target=args.target

    pre_dir_input='../output/method_outputs/'
    raw_data_path='../data/process_data/h5ad_files/'
    save_dir_output='../output/evaluation/'
    
    pre_dir=pre_dir_input+method
    save_dir=save_dir_output+method+'/'
    file_list = os.listdir(pre_dir)

    # read the integrated h5ad file
    file=pre_dir+'/'+target+'.h5ad'
    adata=sc.read_h5ad(file)

    # read the raw data    
    all_target_files = [file for file in file_list if file.endswith('.h5ad')]
    all_raw_files=glob.glob(raw_data_path+target+'_'+'*.h5ad',recursive=True)
    all_raw_files=sorted(all_raw_files)
    print(all_raw_files)
    adata_a=sc.read_h5ad(all_raw_files[0])
    adata_b=sc.read_h5ad(all_raw_files[1])

    # since adata_a, adata_b have already been filtered, normalized and log
    # here we just need to find the hvg and computed trajectory before integration
    sc.pp.highly_variable_genes(adata_a)
    sc.tl.pca(adata_a)
    sc.pp.neighbors(adata_a)
    adata_a.uns["iroot"] = 0
    sc.tl.diffmap(adata_a)
    sc.tl.dpt(adata_a)

    sc.pp.highly_variable_genes(adata_b)
    sc.tl.pca(adata_b)
    sc.pp.neighbors(adata_b)
    adata_b.uns["iroot"] = 0
    sc.tl.diffmap(adata_b)
    sc.tl.dpt(adata_b)

    if method in ['saturn'] and len(all_raw_files)==2:
        # saturn doesn't sort the species' cell name, here we add them back 
        # note: saturn integrates species by the character order
        new_names=adata_a.obs_names.to_list()+adata_b.obs_names.to_list()
        old_names=list(adata.obs_names)
        rename_dict = dict(zip(old_names, new_names))
        adata.obs.rename(index=rename_dict,inplace=True)

    # below is to read the one2one orthologous mapping file to subset the datasets and then merge them
    if target in ['task4']:
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
    elif target in ['task6']:
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

        mapping=mapping.set_index('Stylophora')
        mapping=mapping.loc[adata_b.var_names]
        adata_b.var_names=mapping['gene_name']

        common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))  
        adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
        adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]
        
    elif target in ['task8']:
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

    elif target in ['task9','task9-1']:
        mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/ant_to_mouse.csv')
        mapping=mapping.iloc[:,3:5]
        mapping.drop_duplicates(inplace=True,ignore_index=True)
  
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

    elif target in ['task10']:
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

    elif target in ['task11']:
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

    elif target in ['task12']:
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

    elif target in ['task13']+[f'task{i}' for i in range(33,40)]:
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
        mapping=mapping.loc[mapping['sea_urchin_gene'].isin(adata_a.var_names)] 
        adata_a=adata_a[:,adata_a.var_names.isin(mapping['sea_urchin_gene'])] 
        if method not in ['seurat4','harmony','fastmnn']:
            mapping=mapping.set_index('fish_gene')
            mapping=mapping.loc[adata_b.var_names]
            target_orth_names=mapping['sea_urchin_gene']
            adata_b.var_names=target_orth_names  
        else:
            mapping=mapping.set_index('sea_urchin_gene')
            mapping=mapping.loc[adata_a.var_names]
            target_orth_names=mapping['fish_gene']
            adata_a.var_names=target_orth_names
        common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
        adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
        adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]
    elif target in ['task14']:
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

    elif target in ['task15']:
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

    elif target in ['task16']:
        mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/schmid_to_human.csv')
        mapping.drop_duplicates(inplace=True,ignore_index=True)    
        mapping = mapping.drop_duplicates(subset="human_gene")
        a_overlap=list(set(adata_a.var_names) & set(mapping['human_gene']))
        adata_a=adata_a[:,adata_a.var_names.isin(mapping['human_gene'])]
        mapping=mapping.loc[mapping['human_gene'].isin(a_overlap)] 
        
        b_overlap=list(set(adata_b.var_names) & set(mapping['schmidtea_gene']))
        adata_b=adata_b[:,adata_b.var_names.isin(mapping['schmidtea_gene'])]  ß
        mapping=mapping.loc[mapping['schmidtea_gene'].isin(b_overlap)] 
        mapping=mapping.set_index('schmidtea_gene')
        mapping=mapping.loc[adata_b.var_names]
        target_orth_names=mapping['human_gene']

        adata_b.var_names=target_orth_names   # convert both to the fish gene symbol
        common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
        adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
        adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]
    
    elif target in ["task" + str(i)  for i in range(17,33)]:
        if target in ['task17']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/cat_to_tiger.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            a_name='cat_gene'
            b_name='tiger_gene'
        elif target in ['task18']:
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/cat_to_dog.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            a_name='cat_gene'
            b_name='dog_gene'
        elif target == 'task19':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_to_MF.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            a_name='MF_gene'
            b_name='human_gene'
        elif target == 'task20':
            mapping=pd.read_csv('/ibex/scratch/projects/c2101/benchmark/one2one_orthologs/human_to_MM.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            a_name='MM_gene'
            b_name='human_gene'
        elif target == 'task21':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_to_mouse_mm10.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)    
            a_name='human_gene'
            b_name='mouse_gene'
        elif target == 'task22':
            mapping=pd.read_csv('/ibex/scratch/projects/c2101/benchmark/one2one_orthologs/mouse_to_MF.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='MF_gene'
            b_name='mouse_gene'
        elif target == 'task23':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/mouse_to_MM.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='MM_gene'
            b_name='mouse_gene'
        elif target == 'task24':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_to_pig.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='human_gene'
            b_name='pig_gene'
        elif target == 'task25':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/mouse_to_pig.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='mouse_gene'
            b_name='pig_gene'
        elif target == 'task26':
            mapping=pd.read_csv('/ibex/scratch/projects/c2101/benchmark/one2one_orthologs/pig_to_MF.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='MF_gene'
            b_name='pig_gene'
        elif target == 'task27':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/pig_to_MM.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='MM_gene'
            b_name='pig_gene'
        elif target == 'task28':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/human_to_fish.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='fish_gene'
            b_name='human_gene'
        elif target == 'task29':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/fish_to_MF.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='MF_gene'
            b_name='fish_gene'
        elif target == 'task30':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/fish_to_cat.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='cat_gene'
            b_name='fish_gene'
        elif target == 'task31':
            mapping=pd.read_csv('../OrthoFinder/one2one_orthologs/pig_to_fish.csv')
            mapping.drop_duplicates(inplace=True,ignore_index=True)
            a_name='fish_gene'
            b_name='pig_gene'
        elif target == 'task32':
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
        adata_a=adata_a[:,adata_a.var_names.isin(mapping[a_name])]
        if method in ['seurat4','fastmnn','harmony'] and target in ['task19','task20','task29','task32']:
            mapping = mapping.drop_duplicates(subset=a_name)
            mapping=mapping.set_index(a_name)
            mapping=mapping.loc[adata_a.var_names]
            target_orth_names=mapping[b_name]
            adata_a.var_names=target_orth_names
        else:
            mapping = mapping.drop_duplicates(subset=b_name)
            mapping=mapping.set_index(b_name)
            mapping=mapping.loc[adata_b.var_names]
            target_orth_names=mapping[a_name]
            adata_b.var_names=target_orth_names   
        common_genes=list(set(adata_b.var_names) & set(adata_a.var_names))   
        adata_a=adata_a[:,adata_a.var_names.isin(common_genes)]
        adata_b=adata_b[:,adata_b.var_names.isin(common_genes)]
       
    adata_ls = ad.concat([adata_a,adata_b],join='inner')    

    # double check the adata_ls and adata have same cells
    common_cells=list(set(adata.obs_names) & set(adata_ls.obs_names))
    adata=adata[adata.obs_names.isin(common_cells),:]
    adata_ls=adata_ls[adata_ls.obs_names.isin(common_cells),:]

    if method in ['saturn']:
        if target == 'task4':
            adata.obs['batch']=adata.obs['species_tissue'] # for task4 to check the nested batch
        else:
            adata.obs['batch']=adata.obs['species']
        adata.obs['celltype']=adata.obs['labels2']
        adata.obsm['X_saturn']=adata.X

    if method in ['samap']:
        adata.obs['celltype']=adata.obs['celltype'].apply(lambda x: x[3:])
        adata.obs['batch']=adata.obs['batch']  
        adata.obs['batch']=adata.obs['batch'].apply(lambda x: x[3:])
        all_var_names=pd.DataFrame(adata.var_names)
        all_var_names[0]=all_var_names[0].apply(lambda x: x[3:])
        adata.var_names=all_var_names[0]
    
    adata.obs_names_make_unique()
    adata_ls.obs_names_make_unique()
    adata.obs['batch']=adata.obs['batch'].astype('category')
    adata.obs['celltype']=adata.obs['celltype'].astype('category')
    adata_ls.obs['batch']=adata_ls.obs['batch'].astype('category')
    adata_ls.obs['celltype']=adata_ls.obs['celltype'].astype('category')

    if method in ['scanorama','seurat4']:
        emb_label='X_pca'   # the pca is already computed by the reconstructed data
    elif method in ['fastmnn']:
        emb_label='X_mnn'
    elif method in ['scVI']:
        emb_label='X_scVI'
    elif method in ['scGen']:
        emb_label='corrected_latent'
    elif method in ['saturn']:
        emb_label='X_saturn'
    elif method in ['harmony']:
        emb_label='X_harmony'

    if method not in ['bbknn','samap']:
        sc.pp.neighbors(adata,use_rep=emb_label)

    if method in ['seurat4','scanorama','fastmnn','scGen']:
        ## calculate the hvg consrevation for seurat, scanorama, scGen and fastmnn
        hvg=scib.metrics.hvg_overlap(adata_ls,adata,batch_key='batch')
    else:
        hvg=None

    ## calculate the trajectory conservation 
    traj=scib.me.trajectory_conservation(adata_ls, adata, label_key="celltype")

    # calculate the graph connectivity
    graph_conn=scib.me.graph_connectivity(adata, label_key="celltype")

    # calculate the kbet
    if method in ['seurat4','scanorama']:
        type_label='full'
        kbet=scib.me.kBET(adata, batch_key="batch", label_key="celltype", type_=type_label,embed=emb_label)
    elif method in ['samap','bbknn']:
        type_label='knn'
        kbet=scib.me.kBET(adata, batch_key="batch", label_key="celltype", type_=type_label)
    else:
        type_label='embed'
        kbet=scib.me.kBET(adata, batch_key="batch", label_key="celltype", type_=type_label,embed=emb_label)

    # calculate the PCR batch
    if method in ['seurat4','scanorama']:
        pcr=scib.me.pcr_comparison(adata_ls, adata, covariate="batch")
    elif method not in ['samap','bbknn']:
        pcr=scib.me.pcr_comparison(adata_ls, adata, covariate="batch",embed=emb_label)
    else:
        pcr=None
            
    results={'NMI_cluster/label':None,'ARI_cluster/label':None,'ASW_label':None,'ASW_label/batch':None,'PCR_batch':None,'cell_cycle_conservation':None,'isolated_label_F1':None,'isolated_label_silhouette':None,'graph_conn':None,'kBET':None,'iLISI':None,'cLISI':None,'hvg_overlap':None,'trajectory':None,'true_trajectory':None}

    results['graph_conn']=graph_conn
    results['kBET']=kbet
    results['hvg_overlap']=hvg
    results['trajectory']=traj
    results['PCR_batch']=pcr

    df = pd.DataFrame(list(results.items()), columns=['Variable', 0]).set_index('Variable').rename_axis(None)
    df.to_csv(save_dir+target+'_'+method+'_scib_output.csv')

    
