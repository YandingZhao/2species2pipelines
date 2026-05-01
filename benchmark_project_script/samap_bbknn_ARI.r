# below code is to run the ARI for knn output methods (samap and bbknn)
# to run the code, use the following command:
# Rscript samap_bbknn_ARI.r <target> <method>
library(data.table)
library(mclust)
library('NbClust')

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
        df=read.csv(paste0(pre_dir,target,'_',method,'_leiden_clusters.csv'),header=T)
        rownames(df)<-make.unique(df$X)
        celltype_label='celltype'
        batch_label='batch'

        # run function 20 times, each time extract 80% of data
        percent_extract <- 0.8
        nbiters <- 20
        it <- c()
        total_ari_batch <- c()
        total_ari_celltype <- c()
        myData=df

        nbct=length(unique(df[,celltype_label]))
        ce_types<-unique(myData[,celltype_label])
        # start loop for 20 times
        for(i in 1:nbiters) {
        
            # select cells for the subsampled dataset
            selectedcells<-vector()
            for (g in 1:nbct){
                cellpool<-which(myData[,celltype_label]==ce_types[g])
                ori_nbcells<-length(cellpool)
                cells_extract<-sample(cellpool, size=round(ori_nbcells*percent_extract), replace = F)
                selectedcells<-c(selectedcells, cells_extract)
            }
        
            selectedcells<-sort(selectedcells)
            
            # create the subsampled dataset
            myPCAExt <- myData[selectedcells,]

            # Following clustering, get list of common cell types
            mySample <- myPCAExt[,c(batch_label,celltype_label)]
            print(unique(mySample[,celltype_label]))
            batches <- unique(mySample[,batch_label])
            print(batches)
        
            ctls <- list()
            count <- 0
            for (b in batches){
                count <- count + 1
                ct <- unique(mySample[which(mySample[,batch_label]==b), celltype_label])
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
        
            cells_common <- rownames(mySample)[which(mySample[,celltype_label] %in% ct_common)]
            print(paste("Number of common cells:", length(cells_common)), quote = F)
            
            # create dataset with only common cells
            smallData <- myData[cells_common,]

            # run ARI
            #ari_batch <- mclust::adjustedRandIndex(smallData[,batch_label], smallData$leiden)
            ari_batch_per_celltype <- sapply(ct_common, function(ct) {
                cells_ct <- smallData[smallData[, celltype_label] == ct, ]
                if (nrow(cells_ct) > 1) {
                    mclust::adjustedRandIndex(cells_ct[, batch_label], cells_ct$leiden)
                } else {
                    NA  # Skip if not enough cells for a meaningful ARI
                }
            })
            ari_batch_per_celltype <- ari_batch_per_celltype[!is.na(ari_batch_per_celltype)]
            ari_batch <- mean(ari_batch_per_celltype)
            ari_celltype<-mclust::adjustedRandIndex(myPCAExt[,celltype_label], myPCAExt$leiden)
            it <- c(it,i)
            total_ari_batch <- c(total_ari_batch, ari_batch)
            total_ari_celltype <- c(total_ari_celltype, ari_celltype)
        }

        # once looped 20 times to produce a total of 40 scores, next step follows:
        it <- c(it,nbiters+1)
        total_ari_batch <- c(total_ari_batch, mean(total_ari_batch))
        total_ari_celltype <- c(total_ari_celltype, mean(total_ari_celltype))
        
        methods <- rep(method, nbiters)
        methods <- c(methods,paste0(method,'_mean'))

        # create final dataframe containing raw and median ARI scores
        myARI <- data.frame("use_case"=methods, 
                            "iteration"=it,
                            "ari_batch"=total_ari_batch, 
                            "ari_celltype"=total_ari_celltype)
        
        write.table(myARI,file=paste0(save_dir,target,'_',method,'_ARI.txt'),sep='\t',quote=F,row.names=F,col.names=T)
    }
}
