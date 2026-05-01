#' fastMNN pipeline
#' make sure the batch label and cell type label between two datasets are the same
#' @param merge_rds: the list of the rds object of the two datasets
#' @param batch_label: labels in the metadata for batch in two datasets
#' @param celltype_label: label in the metadata for celltype in two datasets
#' @param save_prefix: the path to the folder that save the output files
#' @param task_name: the name of the task without "_", eg 'task4'
#' @param random_seed: the random seed

# here we run the fastMNN within the seurat v4 workflow
# the output is the corrected embedding and the integrated seurat object which has the reconstructed expression matrix

library(Seurat)
library(SeuratWrappers)
random_seed1 <- sample(1:1000000, 1)

mnn_preprocess <- function(merge_rds,batch_label='batch',celltype_label='celltype',
                                save_prefix='../output/method_outputs/fastMNN/',task_name,random_seed=random_seed1){
    #random_seed <- sample(1:1000000, 1)
    print(random_seed)
    set.seed(random_seed)
    start_time=proc.time()    
    ## create seurat object
    b_seurat<-merge_rds
    b_seurat=NormalizeData(b_seurat)
    nfeature_max=dim(b_seurat)[1]
    if (nfeature_max > 2000) {
        nfeature=2000
    } else {
        nfeature=nfeature_max
    }
    b_seurat<-FindVariableFeatures(b_seurat,selection.method = "vst", nfeatures = nfeature)
    b_seurat<-ScaleData(b_seurat)
  
    ## plotting
    k_seed=10
    umap_perplex=30
    umapplot_filename='_fastmnn_umap'
    obj_filename='_fastmnn_obj'
    pca_filename='_fastmnn_pca'
    set.seed(random_seed)
    b_seurat=RunFastMNN(SplitObject(b_seurat,split.by=batch_label))
    
    fastmnn_res <- as.data.frame(b_seurat@reductions$mnn@cell.embeddings)
    cells_use <- rownames(fastmnn_res)
    fastmnn_res$batchlb <- b_seurat@meta.data[, batch_label]
    fastmnn_res$celltype <- b_seurat@meta.data[, celltype_label]
    end_time=proc.time()
    elapsed_time <- end_time - start_time
    log_file <- "../output/method_outputs/fastMNN/fastMNN_timing_log.txt"
    write(paste("Task name:", task_name, "Seed:", random_seed, "User time:", elapsed_time["user.self"], "System time:", elapsed_time["sys.self"], "Elapsed time:", elapsed_time["elapsed"], "\n"), file = log_file, append = TRUE)
    ## save the result
    write.table(fastmnn_res, file=paste0(save_prefix,task_name, "_",random_seed,"_fastmnn_embedding.txt"), quote=F, sep='\t', row.names = T, col.names = NA)
    print('embedding table saved!')
    saveRDS(b_seurat, file=paste0(save_prefix,task_name,"_",random_seed,"_fastmnn_integration.rds"))
    print('integrated rds file saved!')

    ## save the UMAP plot
    b_seurat <- RunUMAP(b_seurat, reduction = "mnn", dims = 1:30)
    p11 <- DimPlot(object = b_seurat, reduction = 'umap', group.by = batch_label)
    p12 <- DimPlot(object = b_seurat, reduction = 'umap', group.by = celltype_label)

    library(cowplot)
    png(paste0(save_prefix,task_name,"_fastmnn_umap_plot.png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p11, p12))
    dev.off()

    pdf(paste0(save_prefix,task_name,"_fastmnn_umap_plot.pdf"),width=15,height=7,paper='special') 
    print(plot_grid(p11, p12))
    dev.off()
    return(b_seurat)
    print('umap plot saved!')
}
