# below code is to run the LISI metric for the knn output methods (samap and bbknn)
# to run the code, use the following command:
# python samap_bbknn_evaluate_LISI.py --all_targets target1 target2 target3 --methods samap bbknn
import pandas as pd
import numpy as np
import scanpy as sc
import os
import sys
import scipy.sparse
import argparse
sys.path.insert(1, './core_script')
os.chdir('./core_script')
from evaluation_LISI_knn import *

# read the knn method(bbknn,samap) output data 
pre_dir_input='../../output/method_outputs/'
save_dir_output='../../output/evaluation/'

def main():
    parser = argparse.ArgumentParser(description='samap and bbknn output evaluation NMI')
    parser.add_argument('--all_targets', nargs='+', help='list of target task names')
    parser.add_argument('--methods', nargs='+', help='list of methods')
    args = parser.parse_args()
    methods = args.methods
    all_targets = args.all_targets
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
            # run lisi
            lisi_output=evaluate_LISI_knn(adata,batch_key,label_key)
            means=lisi_output.mean(axis=0)
            lisi_output.loc['mean']=means
            lisi_output.to_csv(save_dir+target+'_'+method+'_LISI.csv')

main()

        
