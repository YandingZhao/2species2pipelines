#' Hrmony(0.1.1) pipeline
#' make sure the batch label and cell type label between two datasets are the same
#' @param merge_rds: list of the targets datasets rds object
#' @param batch_label: labels in the metadata for batch in two datasets
#' @param celltype_label: label in the metadata for celltype in two datasets
#' @param save_prefix: the path to the folder that save the output files
#' @param task_name: the name of the task without "_", eg 'task4'
#' @param random_seed: the random seed
# here we run harmony within the seurat v4 workflow with the maximal number of clusters (50) and the maximal number of iterations (100).
# the output is the corrected PCA space

library(Seurat)
library(harmony)
random_seed1 <- sample(1:1000000, 1)

harmony_preprocess <- function(merge_rds,batch_label='batch',celltype_label='celltype',
                                save_prefix='../output/method_outputs/harmony/',task_name,random_seed=random_seed1){
    
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

    ## harmony setting
    theta_harmony=2
    numcluster=50
    max_iter_cluster=100

    ## plotting
    k_seed=10
    umap_perplex=30
    umapplot_filename='_harmony_umap'
    obj_filename='_harmony_obj'
    pca_filename='_harmony_pca'

    b_seurat <- RunPCA(b_seurat, npcs = 30, verbose = FALSE)
    b_seurat <- RunHarmony(b_seurat, batch_label,theta=theta_harmony,plot_convergence=T,
                            nclust=numcluster,max.iter.harmony=max_iter_cluster)
    
    harmony_res <- as.data.frame(b_seurat@reductions$harmony@cell.embeddings)
    cells_use <- rownames(harmony_res)
    harmony_res$batchlb <- b_seurat@meta.data[, batch_label]
    harmony_res$celltype <- b_seurat@meta.data[, celltype_label]
    
    end_time=proc.time()
    elapsed_time <- end_time - start_time
    log_file <- "../output/method_outputs/harmony/harmony_timing_log.txt"
    write(paste("Task name:", task_name, "Seed:", random_seed, "User time:", elapsed_time["user.self"], "System time:", elapsed_time["sys.self"], "Elapsed time:", elapsed_time["elapsed"], "\n"), file = log_file, append = TRUE)

    ## save the result
    write.table(harmony_res, file=paste0(save_prefix,task_name, "_",random_seed,"_harmony_pca.txt"), quote=F, sep='\t', row.names = T, col.names = NA)
    print('pca table saved!')
    saveRDS(b_seurat, file=paste0(save_prefix,task_name, "_",random_seed,"_harmony_integration.rds"))
    print('integrated rds file saved!')

    ## save the UMAP plot
    b_seurat <- RunUMAP(b_seurat, reduction = "harmony", dims = 1:30)
    p11 <- DimPlot(object = b_seurat, reduction = 'umap', group.by = batch_label)
    p12 <- DimPlot(object = b_seurat, reduction = 'umap', group.by = celltype_label)

    library(cowplot)
    png(paste0(save_prefix,task_name,"_harmony_umap_plot.png"),width = 2*1000, height = 800, res = 2*72)
    print(plot_grid(p11, p12))
    dev.off()

    pdf(paste0(save_prefix,task_name,"_harmony_umap_plot.pdf"),width=15,height=7,paper='special') 
    print(plot_grid(p11, p12))
    dev.off()
    return(b_seurat)
    print('umap plot saved!')
}
