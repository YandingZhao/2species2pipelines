# this script is to run the ASW and NMI metrics for the embedding/reconstruction output methods
# to run the script, use the following command:
# python harmony_ASW_NMI.py --all_targets target1 target2 target3 --methods method1 method2 method3

import os
import sys
sys.path.insert(1, './core_script')
from evaluation_ASW_NMI import *
import argparse
os.chdir('./core_script')

pre_dir_input='../../output/method_outputs/'
save_dir_output='../../output/evaluation/'

def main():
    parser = argparse.ArgumentParser(description='evaluation for ASW and NMI for all tasks')
    parser.add_argument('--all_targets', nargs='+', help='list of target task names')
    parser.add_argument('--methods', nargs='+', help='list of methods')
    args = parser.parse_args()
    all_targets = args.all_targets
    methods = args.methods
    for method in methods:
        print(method) 
        for target in all_targets:
            #print(target)
            if method in ['harmony','scanorama','seurat4']:
                file_name=target+'_'+method+"_pca.txt"
            elif method in ['scVI','fastmnn','scGen','saturn','saturn_one2one']:
                file_name=target+'_'+method+'_embedding.txt'
            pre_dir=pre_dir_input+method
            save_dir=save_dir_output+method+'/'
            if method in ['scanorama','scVI','scGen','saturn','saturn_one2one']:
                sep_symbol=','
            else:
                sep_symbol='\t'
            adata=createAnnData(pre_dir,file_name,sep_symbol)
            task_name=target
            print(task_name)
            
            adata.obs['batch']=adata.obs['batch'].astype(str)
            AWS_df=silhouette_coeff_ASW(adata,method_use=method,save_dir=save_dir, task_name=task_name, percent_extract=0.8)
            nmi_dict={}
            nmi_dict['nmi_value_celltype']=nmi(adata,'cell_type')
            nmi_dict['nmi_value_batch']=nmi(adata,'batch')
            nmi_df=pd.DataFrame.from_dict([nmi_dict])
            nmi_df.to_csv(save_dir+task_name+'_'+method+'_NMI.csv')

main()