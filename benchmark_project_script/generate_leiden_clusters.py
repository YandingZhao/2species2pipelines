## this code is to generate the best clustering results of the leiden algorithm for bbknn and samap output based on the ari (celltype) score
## the output is saved in the output/method_outputs/{method}/{target}_{method}_leiden_clusters.csv
# to run this code, you need to run the following command:
# python generate_leiden_clusters.py --all_targets task3 task4 --methods bbknn samap
import pandas as pd
import numpy as np
import scanpy as sc
import os
import sys
import scipy.sparse
from sklearn.metrics.cluster import adjusted_rand_score
import argparse


def optimal_leiden(adata,label_key):
    '''
    label_key: celltype or cell_lineage
    adata: adata object contain the graph integration output
    '''
    n=20
    resolution = [2*x/n for x in range(1,n+1)]
    best_ari=0
    for res in resolution:
        sc.tl.leiden(adata,resolution=res)
        ari=adjusted_rand_score(adata.obs[label_key],adata.obs['leiden'])
        if ari>best_ari:
            best_ari=ari
            adata.obs['leiden_best']=adata.obs['leiden']
        del adata.obs['leiden']
    return adata


def main():
    parser = argparse.ArgumentParser(description='run the leiden clustering algorithm on samap and bbknn')
    parser.add_argument('--all_targets', nargs='+', help='list of target task names')
    parser.add_argument('--methods', nargs='+', help='list of methods')
    args = parser.parse_args()
    all_targets = args.all_targets
    methods = args.methods

    pre_dir_input='../output/method_outputs/'
    save_dir_output='../output/method_outputs/'

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

            # run leiden
            optimal_leiden(adata,label_key)
            df=adata.obs[[batch_key,label_key,'leiden_best']]
            df=adata.obs[[batch_key,label_key,'leiden']]
            # rename the column 'leiden_best' to 'leiden'
            df.rename(columns={'leiden_best':'leiden'},inplace=True)
            df.to_csv(save_dir+target+'_'+method+'_leiden_clusters.csv')

main()