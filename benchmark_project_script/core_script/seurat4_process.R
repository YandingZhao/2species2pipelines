#' Seurat(4.3.0) pipeline
#' make sure the batch label and cell type label between two datasets are the same
#' @param target_list: input data list including a and b datasets
#' @param batch_label: labels in the metadata for batch in two datasets
#' @param celltype_label: label in the metadata for celltype in two datasets
#' @param save_prefix: the path to the folder that save the output files
#' @param task_name: the name of the task without "_", eg 'task4'
#' @param random_seed: the random seed
# the output would be the 30 PCs 

library(Seurat)
random_seed1 <- sample(1:1000000, 1)

seurat4_preprocess <- function(target_list,batch_label='batch',celltype_label='celltype',save_prefix='../output/method_outputs/seurat4/',task_name,random_seed=random_seed1){

    print(random_seed)
    set.seed(random_seed)
    start_time=proc.time()
    ## perform the integration
    #target_list=list(target_a,target_b)
    target_list <- lapply(X = target_list, FUN = function(x) {
        x <- NormalizeData(x)
        nfeature_max=min(dim(target_list[[1]])[1],dim(target_list[[2]])[1])
        if (nfeature_max > 2000) {
            nfeature=2000
        } else {
            nfeature=nfeature_max
        }
        # here we set the number of features to be selected as 2000 if the number of features in the dataset is larger than 2000, else is the number of features in the dataset
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeature)
    })
    features <- SelectIntegrationFeatures(object.list = target_list)
    target_anchors <- FindIntegrationAnchors(object.list = target_list, anchor.features = features)

    target_combined <- IntegrateData(anchorset = target_anchors)
    DefaultAssay(target_combined) <- "integrated"
    target_combined <- ScaleData(target_combined, verbose = FALSE)
    target_combined <- RunPCA(target_combined, npcs = 30, verbose = FALSE)
    target_combined <- RunUMAP(target_combined, reduction = "pca", dims = 1:30)
    
    seurat4_res <- as.data.frame(target_combined@reductions$pca@cell.embeddings)
    cells_use <- rownames(seurat4_res)
    seurat4_res$batchlb <- target_combined@meta.data[, batch_label]
    seurat4_res$celltype <- target_combined@meta.data[, celltype_label]

    end_time=proc.time()
    elapsed_time <- end_time - start_time
    log_file <- "../output/method_outputs/seurat4/seurat4_timing_log.txt"
    write(paste("Task name:", task_name, "Seed:", random_seed, "User time:", elapsed_time["user.self"], "System time:", elapsed_time["sys.self"], "Elapsed time:", elapsed_time["elapsed"], "\n"), file = log_file, append = TRUE)

    ## save the result
    write.table(seurat4_res, file=paste0(save_prefix,task_name,"_",random_seed, "_seurat4_pca.txt"), quote=F, sep='\t', row.names = T, col.names = NA)
    print('pca table saved!')
    saveRDS(target_combined, file=paste0(save_prefix,task_name,"_",random_seed,"_seurat4_integration.rds"))
    print('integrated rds file saved!')

    ## save the UMAP plot
    p11 <- DimPlot(object = target_combined, reduction = 'umap', group.by = batch_label)
    p12 <- DimPlot(object = target_combined, reduction = 'umap', group.by = celltype_label)

    library(cowplot)
    png(paste0(save_prefix,task_name,"_seurat4_umap_plot.png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p11, p12))
    dev.off()

    pdf(paste0(save_prefix,task_name,"_seurat4_umap_plot.pdf"),width=15,height=7,paper='special') 
    print(plot_grid(p11, p12))
    dev.off()
    return(target_combined)
    print('umap plot saved!')
}
