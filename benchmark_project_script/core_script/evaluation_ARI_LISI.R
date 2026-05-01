### this code is to run the ARI and LISI evaluation for the method outputs

library(data.table)
library(NbClust)
library(mclust)

ari_calcul_sampled <- function(myData, cpcs, isOptimal=FALSE, 
                               method_use='resnet',
                               base_name='', maxiter=30, 
					           celltypelb='celltype', batchlb='batchlb' )
{

    set.seed(0)
    
    # get number of unique cell types
    nbct <- length(unique(myData[,celltypelb]))
    
    # get vector of unique cell types
    ce_types<-unique(myData[,celltypelb])
    
    # run function 20 times, each time extract 80% of data
    nbiters <- 20
    percent_extract <- 0.8
    
    it <- c()
    total_ari_batch <- c()
    total_ari_celltype <- c()
    
    # start loop for 20 times
    for(i in 1:nbiters) {
        
        # select cells for the subsampled dataset
        selectedcells<-vector()
        for (g in 1:nbct){
          cellpool<-which(myData[,celltypelb]==ce_types[g])
          ori_nbcells<-length(cellpool)
          cells_extract<-sample(cellpool, size=round(ori_nbcells*percent_extract), replace = F)
          selectedcells<-c(selectedcells, cells_extract)
        }
    
        selectedcells<-sort(selectedcells)
        
        # create the subsampled dataset
        myPCAExt <- myData[selectedcells,]
        
        ###############################
        # Clustering
        ###############################
        
        if(!isOptimal){  # isOptimal==FALSE
          # nbct : k equal number of unique cell types in the dataset
          clustering_result <- kmeans(x = myPCAExt[,cpcs], centers=nbct, iter.max = maxiter)
          myPCAExt$clusterlb <- clustering_result$cluster
        
        } else if(isOptimal){
          nbclust_result<-NbClust(data=myPCAExt[,cpcs], method = "kmeans", 
                                  min.nc = max(nbct-2, 2), max.nc = nbct+4)
          myPCAExt$clusterlb <-nbclust_result$Best.partition	  
        }
        
        # assign the current myPCAExt to a unique object so that it can be stored later 
        assign(paste0("myPCAExt",i), myPCAExt)
        print('Nb clusters: ')
        print(length(unique(myPCAExt$clusterlb)))
        
        
        # Following clustering, get list of common cell types
        mySample <- subset(myPCAExt,select=c('celltype', 'batchlb'))
        print(unique(mySample[,celltypelb]))
        batches <- unique(mySample[,batchlb])
        print(batches)
        
        ctls <- list()
        count <- 0
        for (b in batches){
        count <- count + 1
        ct <- unique(mySample[which(mySample[,batchlb]==b), celltypelb])
        ctls[[count]] <- ct
        }
        
        for(t in rep(1:length(ctls))){
        if(t==1){
            ct_common <- intersect(ctls[[t]], ctls[[t+1]])    
        }
        if(t>2){   #more than 2 batches
            ct_common <- intersect(ct_common, ctls[[t]])    
        }
        }
        # ct_common: common cell types amongst all batches
        
        cells_common <- rownames(mySample)[which(mySample[,celltypelb] %in% ct_common)]
        print(paste("Number of common cells:", length(cells_common)), quote = F)
        
        # create dataset with only common cells
        smallData <- myPCAExt[cells_common,]
        #assign(paste0("smallData",i), smallData)
        
        ###############################
        # ARI
        ###############################
        
        # run ARI
        #ari_batch <- mclust::adjustedRandIndex(smallData[,batchlb], smallData$clusterlb)
        # Calculate ARI for each cell type and average them for ARI batch
        ari_batch_per_celltype <- sapply(ct_common, function(ct) {
            cells_ct <- smallData[smallData[, celltypelb] == ct, ]
            if (nrow(cells_ct) > 1) {
                mclust::adjustedRandIndex(cells_ct[, batchlb], cells_ct$clusterlb)
            } else {
                NA  # Skip if not enough cells for a meaningful ARI
            }
        })
        ari_batch_per_celltype <- ari_batch_per_celltype[!is.na(ari_batch_per_celltype)]
        ari_batch <- mean(ari_batch_per_celltype)

        ari_celltype<-mclust::adjustedRandIndex(myPCAExt[,celltypelb], myPCAExt$clusterlb)
        
        it <- c(it,i)
        total_ari_batch <- c(total_ari_batch, ari_batch)
        total_ari_celltype <- c(total_ari_celltype, ari_celltype)
    }   # End of loop 
    
    
    # once looped 20 times to produce a total of 40 scores, next step follows:
    it <- c(it,nbiters+1)
    total_ari_batch <- c(total_ari_batch, mean(total_ari_batch))
    total_ari_celltype <- c(total_ari_celltype, mean(total_ari_celltype))
    
    methods <- rep(method_use, nbiters)
    methods <- c(methods,paste0(method_use,'_mean'))
    
    # create final dataframe containing raw and median ARI scores
    myARI <- data.frame("use_case"=methods, 
                        "iteration"=it,
                        "ari_batch"=total_ari_batch, 
                        "ari_celltype"=total_ari_celltype)
    
    # write final dataframe to a text file
    write.table(myARI, file = paste0(base_name,method_use,"_ARI.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep="\t")
    
    print('Save output in folder')
    print(base_name)
    
    return(list(myARI, myPCAExt1, myPCAExt2, myPCAExt3, myPCAExt4, myPCAExt5, myPCAExt6,
                myPCAExt7, myPCAExt8, myPCAExt9, myPCAExt10, myPCAExt11, 
                myPCAExt12,  myPCAExt13, myPCAExt14, myPCAExt15, myPCAExt16,  
                myPCAExt17, myPCAExt18, myPCAExt19, myPCAExt20))
    }


##########################################################################
## calculate the LISI score
get_celltype_common <- function(myData){
    print(dim(myData))
    batches <- unique(myData$batchlb)
    celltypels <- unique(myData$celltype)
    print(celltypels)
    if(length(batches)>1){
        ctls <- list()
        count <- 0
        for (b in batches){
        count <- count + 1
        ct <- unique(myData[which(myData$batchlb==b),'celltype'])
        ctls[[count]] <- ct
        print("Batch")
        print(b)
        print(ct)
        }
    for(i in rep(1:length(ctls))){
      if(i==2){
        ct_common <- intersect(ctls[[i-1]], ctls[[i]])    
      }
      if(i>2){   #more than 2 batches
        ct_common <- intersect(ct_common, ctls[[i]])    
      }
    }
    ct_common <- unique(ct_common)
    print(ct_common)
    cells_common <- rownames(myData)[which(myData$celltype %in% ct_common) ]
    return(list('ct_common'=ct_common, 'cells_common'=cells_common, 'batches'=batches, 'celltypels'=celltypels))
    } else{
        return(NULL)
    } 
    
}


run_LISI_final <- function(myData, plx=40,save_dir,task_name,method_used){
  
    myPCA <- myData

    cpcs <- colnames(myPCA)[1:(ncol(myPCA)-2)]
    lisi_embeddings <- myPCA[,cpcs]
    
    colnames(myPCA)[grep('[cC]ell_?[tT]ype',colnames(myPCA))] <- 'cell_type'
    colnames(myPCA)[grep('([bB]atch)|(BATCH)|(batchlb)',colnames(myPCA))] <- 'batch'
    
    lisi_meta_data <- subset(myPCA, select=c('batch','cell_type'))
    
    lisi_label = c('batch', 'cell_type')
    
    lisi_res <- lisi::compute_lisi(lisi_embeddings, lisi_meta_data, lisi_label,perplexity = plx)
    lisi_res$cell <- rownames(lisi_embeddings)
    
    lisi_batch <- subset(lisi_res,select=c('batch','cell'))
    lisi_celltype <- subset(lisi_res,select=c('cell_type','cell'))
    lisi_batch_mean=mean(lisi_batch$batch)
    lisi_celltype_mean=mean(lisi_celltype$cell_type)
    lisi_batch[nrow(lisi_batch)+1,]=c(lisi_batch_mean,'all_mean')
    lisi_celltype[nrow(lisi_celltype)+1,]=c(lisi_celltype_mean,'all_mean')
    write.table(lisi_batch, paste0(save_dir,task_name,'_',method_used,'_lisi_batch_40.txt'), quote=F, sep='\t', row.names=T, col.names=NA)
    write.table(lisi_celltype,paste0(save_dir,task_name,'_',method_used,'_lisi_celltype_40.txt'),quote=F, sep='\t', row.names=T, col.names=NA)

}


#############-------------------------------------------------------#############
run_all_metric <- function(task_name,method_used,pre_dir='../output/method_outputs/seurat4/',save_dir='../output/evaluation/'){
    if (method_used %in% c('harmony','scanorama','seurat4')){
      if (method_used %in% c('scanorama')) {
        sep_symbol=','
      } else {
        sep_symbol='\t'
      }
      pca_df=read.table(paste0(pre_dir,task_name,"_",method_used,'_pca.txt'),sep=sep_symbol,header=T,row.names=1,quote='')
    } else if (method_used %in% c('scVI','fastmnn','scGen','saturn','saturn_one2one')){
      if (method_used %in% c('scVI','scGen','saturn','saturn_one2one')){
        sep_symbol=','
      } else {
        sep_symbol='\t'
      }
      print(paste0(pre_dir,task_name,"_",method_used,'_embedding.txt'))
      pca_df=read.table(paste0(pre_dir,task_name,"_",method_used,'_embedding.txt'),sep=sep_symbol,header=T,row.names=1,quote='')
    }  
    
    colPCA <- colnames(pca_df)[1:(ncol(pca_df)-2)]
    #colPCA=colnames(pca_df)[1:30]

    base_name=paste0(save_dir,task_name,'_')

    # ari metric
    temp<-ari_calcul_sampled(myData=pca_df, cpcs=colPCA, isOptimal=FALSE, 
                        method_use = method_used,  
                        base_name=base_name,celltypelb='celltype', batchlb='batchlb')
    print('ari metric saved!')

    # lisi metric
    run_LISI_final(pca_df,plx=40,save_dir=save_dir,task_name=task_name,method_used=method_used)
    print('lisi metric saved!')

}