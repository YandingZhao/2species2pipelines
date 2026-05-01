# this script is to run the ARI and LISI evaluation for the embedding/reconstruction output methods (except the knn- graph output methods)
# to run the script, use the following command:
# Rscript emb_recon_ARI_LISI.R <target> <method_name>
# note the target should be "task3_" with "_" 

library(Seurat)
library(SeuratObject)
library(Matrix)
library(data.table)

source('./core_script/evaluation_ARI_LISI.R')

args <- commandArgs(trailingOnly=TRUE)

targets <- c(args[1])
methods=c(args[2])
preflex_input_dir='../output/method_outputs/'
preflex_save_dir='../output/evaluation/'

for (method in methods) {
    print(method)
    for (target in targets) {
        print(target)
        pre_dir=paste0(preflex_input_dir,method,'/')
        save_dir=paste0(preflex_save_dir,method,'/')
        print(pre_dir)
        run_all_metric(target,method,pre_dir=pre_dir,save_dir=save_dir)
    }
}

