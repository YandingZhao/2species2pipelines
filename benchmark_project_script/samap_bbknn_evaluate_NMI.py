# below code is to run the NMI evaluation for samap and bbknn
# to run the code, use the following command:
# python samap_bbknn_evaluate_NMI.py --all_targets target1 target2 target3 --methods samap bbknn

import os
import sys
sys.path.insert(1, './core_script')
os.chdir('./core_script')
from evaluation_ASW_NMI import *
import argparse
import scanpy as sc
#methods=['samap','bbknn']
pre_dir_input='../output/method_outputs/'
save_dir_output='../output/evaluation/'

def main():
    parser = argparse.ArgumentParser(description='samap and bbknn output evaluation on NMI')
    parser.add_argument('--all_targets', nargs='+', help='list of target task names')
    parser.add_argument('--methods', nargs='+', help='list of methods')
    args = parser.parse_args()
    all_targets = args.all_targets
    methods = args.methods
    for method in methods:
        print(method) 
        for target in all_targets:
            print(target)
            pre_dir=pre_dir_input+method+'/'
            save_dir=save_dir_output+method+'/'
            file_name=target+'.h5ad'
            adata=sc.read_h5ad(pre_dir+file_name)
            batch_key='batch'   
            label_key='celltype'
            adata.obs[batch_key]=adata.obs[batch_key].astype('category')
            adata.obs[label_key]=adata.obs[label_key].astype('category')
            nmi_dict={}
            nmi_dict['nmi_value_celltype']=nmi(adata,label_key)
            nmi_dict['nmi_value_batch']=nmi(adata,batch_key)
            nmi_df=pd.DataFrame.from_dict([nmi_dict])
            nmi_df.to_csv(save_dir+target+'_'+method+'_NMI.csv')
main()
