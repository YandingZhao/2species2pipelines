library(Seurat)
library(SeuratObject)
library(Matrix)
library(data.table)

### this script is to run the seurat v4 CCA for the integration
### the usage is to run the script in the command line with the following command
### Rscript seurat4_for_all_tasks.R 'task3_' 

### the script will read the data from the folder ../data/extra_data_more_replicates/raw_data/ which store the RDS file for the raw data
### the output data would be store in the folder ../output/method_outputs/seurat4/
setwd('./script')
source('./core_script/seurat4_process.R')

RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}

options(future.globals.maxSize = 8000 * 1024 ^ 2)

args <- commandArgs(trailingOnly=TRUE)

targets <- c(args[1])

for (target in targets){
    print(target)
    celltype_label='celltype'
    batch_label='batch'

    # read the data    
    all_files_input=list.files('../data/process_data/rds_files/',pattern=target)
    input_a=readRDS(paste0('../data/process_data/rds_files/',all_files_input[1]))
    input_b=readRDS(paste0('../data/process_data/rds_files/',all_files_input[2]))
    nfeatures_a <- Matrix::colSums(x = input_a@assays$RNA@counts > 0)
    num.cells_a <- Matrix::rowSums(x = input_a@assays$RNA@counts > 0)
    nfeatures_b <- Matrix::colSums(x = input_b@assays$RNA@counts > 0)
    num.cells_b <- Matrix::rowSums(x = input_b@assays$RNA@counts > 0)
    input_a_rds <- input_a[which(x=num.cells_a>=10),which(x=nfeatures_a>=200)]
    input_b_rds <- input_b[which(x=num.cells_b>=10),which(x=nfeatures_b>=200)]

    # read the one2one orthologs mapping table and then rename the genes names in input_b_rds to the same as input_a_rds
    if (target %in% c('task4_','task5_','task4-1_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/human_mouse.csv',header=T)
        mapping=mapping[,c(4,5)]
        mapping=mapping[!duplicated(mapping),]    

        # remove all the duplicated genes in mouse and human in the mapping file
        mapping=mapping[!duplicated(mapping$mouse_gene),]
        mapping=mapping[!duplicated(mapping$human_gene),]
        
        ################## convert the mouse gene to human orthologous #############
        human_overlap=intersect(rownames(input_a_rds),mapping$human_gene)
        input_a_rds=input_a_rds[mapping$human_gene,]  
        mapping=mapping[mapping$human_gene %in% human_overlap,]  
        input_b_rds=input_b_rds[mapping$mouse_gene,]  
        mapping=mapping[mapping$mouse_gene %in% rownames(input_b_rds),]

        mapping=mapping[match(rownames(input_b_rds),mapping$mouse_gene),]

        # RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes)
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=mapping$human_gene)
    } else if (target %in% c('task6_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/nema_sty.csv',header=T)
        mapping=mapping[!duplicated(mapping),]
        nema_overlap=intersect(rownames(input_a_rds),mapping$gene_name)
        input_a_rds=input_a_rds[mapping$gene_name,]  
        mapping=mapping[mapping$gene_name %in% nema_overlap,]  
        target_sty_names=gsub('.{2}$','',rownames(input_b_rds))
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_sty_names)  

        input_b_rds=input_b_rds[mapping$Stylophora,] 
        mapping=mapping[mapping$Stylophora %in% rownames(input_b_rds),]  

        # order the mouse gene order in the mapping file
        mapping=mapping[match(rownames(input_b_rds),mapping$Stylophora),]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=mapping$gene_name)

        common_genes=intersect(rownames(input_b_rds),rownames(input_a_rds))   #5919
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]
        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))

    } else if (target %in% c('task8_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/fish_to_frog.csv',header=T)
        mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping$fish_gene),]
        mapping=mapping[!duplicated(mapping$frog_gene),]
        frog_overlap=intersect(rownames(input_b_rds),mapping$frog_gene)
        input_b_rds=input_b_rds[mapping$frog_gene,]  
        mapping=mapping[mapping$frog_gene %in% frog_overlap,]

        input_a_rds=input_a_rds[mapping$fish_gene,]
        #fish_overlap=intersect(rownames(input_a_rds),mapping$fish_gene)
        #mapping=mapping[mapping$fish_gene %in% fish_overlap,]

        mapping=mapping[match(rownames(input_b_rds),mapping$frog_gene),]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=mapping$fish_gene)

        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))

        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))  
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]
    } else if (target %in% c('task9_','task9-1_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/ant_to_mouse.csv',header=T)
        mapping=mapping[,c(4,5)]
        mapping=mapping[!duplicated(mapping),]    

        # remove all the duplicated genes in mouse and human in the mapping file    
        mapping=mapping[!duplicated(mapping$mouse_gene),]
        mapping=mapping[!duplicated(mapping$ant_gene),]

        mapping$ant_gene<-gsub('_','-',mapping$ant_gene)

        ################## convert the ant gene to mouse orthologous #############
        ant_overlap=intersect(rownames(input_a_rds),mapping$ant_gene)
        input_a_rds=subset(input_a_rds,features=ant_overlap)
        mapping=mapping[mapping$ant_gene %in% ant_overlap,]  
        mapping=mapping[match(rownames(input_a_rds),mapping$ant_gene),]

        # RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes)
        input_a_rds=RenameGenesSeurat(input_a_rds,newnames=mapping$mouse_gene)

        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]
        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task10_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/fish_to_fly.csv',header=T)
        mapping=mapping[!duplicated(mapping),]   

        ################## convert the fish gene to fly orthologous #############
        fish_overlap=intersect(rownames(input_a_rds),mapping$fish_gene)
        input_a_rds=input_a_rds[mapping$fish_gene,]  
        mapping=mapping[mapping$fish_gene %in% fish_overlap,]  
        fly_overlap=intersect(rownames(input_b_rds),mapping$fly_symbol)
        input_b_rds=input_b_rds[mapping$fly_symbol,]  
        mapping=mapping[mapping$fly_symbol %in% fly_overlap,] 

        target_orth_fish_names=mapping$fish_gene
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_fish_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))  
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task11_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/trichoplax_to_human.csv',header=T)
        mapping=mapping[!duplicated(mapping),]    
        human_overlap=intersect(rownames(input_a_rds),mapping$human_gene)
        input_a_rds=input_a_rds[mapping$human_gene,]  

        mapping=mapping[mapping$human_gene %in% human_overlap,]  
        trichoplax_overlap=intersect(rownames(input_b_rds),mapping$trichoplax_gene)
        input_b_rds=input_b_rds[mapping$trichoplax_gene,]  
        mapping=mapping[mapping$trichoplax_gene %in% trichoplax_overlap,] 

        target_orth_names=mapping$human_gene
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task12_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/ciona_to_nema.csv',header=T)
        mapping=mapping[!duplicated(mapping),]    
        ciona_overlap=intersect(rownames(input_a_rds),mapping$ciona_gene)
        input_a_rds=input_a_rds[mapping$ciona_gene,]  

        mapping=mapping[mapping$ciona_gene %in% ciona_overlap,]  
        nema_overlap=intersect(rownames(input_b_rds),mapping$nema_gene)
        input_b_rds=input_b_rds[mapping$nema_gene,] 
        mapping=mapping[mapping$nema_gene %in% nema_overlap,] 

        target_orth_names=mapping$ciona_gene
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task13_','task33_','task34_','task35_','task36_','task37_','task38_','task39_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/sea_urchin_to_zebrafish.csv',header=T)
        mapping=mapping[!duplicated(mapping),]    
        b_overlap=intersect(rownames(input_b_rds),mapping$fish_gene)
        input_b_rds=input_b_rds[mapping$fish_gene,]  

        mapping=mapping[mapping$fish_gene %in% b_overlap,]  
        a_overlap=intersect(rownames(input_a_rds),mapping$sea_urchin_gene)
        input_a_rds=input_a_rds[mapping$sea_urchin_gene,]  
        mapping=mapping[mapping$sea_urchin_gene %in% a_overlap,] 

        target_orth_names=mapping$fish_gene
        input_a_rds=RenameGenesSeurat(input_a_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task14_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/celegan_to_human.csv',header=T)
        mapping=mapping[!duplicated(mapping),]    
        a_overlap=intersect(rownames(input_a_rds),mapping$celegan_gene)
        input_a_rds=input_a_rds[mapping$celegan_gene,]  
        mapping=mapping[mapping$celegan_gene %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping$human_gene)
        input_b_rds=input_b_rds[mapping$human_gene,]  
        mapping=mapping[mapping$human_gene %in% b_overlap,] 
        mapping=mapping[!duplicated(mapping[,c('human_gene')]),]
        target_orth_names=mapping$celegan_gene
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))  
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task15_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/octopus_to_human.csv',header=T)
        mapping=mapping[!duplicated(mapping),]    
        a_overlap=intersect(rownames(input_a_rds),mapping$human_gene)
        input_a_rds=input_a_rds[mapping$human_gene,] 
        mapping=mapping[mapping$human_gene %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping$octopus_gene)
        input_b_rds=input_b_rds[mapping$octopus_gene,]  
        mapping=mapping[mapping$octopus_gene %in% b_overlap,]
        mapping=mapping[!duplicated(mapping[,c('human_gene')]),] 
        target_orth_names=mapping$human_gene
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))  
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task16_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/schmid_to_human.csv',header=T)
        mapping=mapping[!duplicated(mapping),]    
        a_overlap=intersect(rownames(input_a_rds),mapping$human_gene)
        input_a_rds=input_a_rds[mapping$human_gene,]  

        mapping=mapping[mapping$human_gene %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping$schmidtea_gene)
        input_b_rds=input_b_rds[mapping$schmidtea_gene,]  
        mapping=mapping[mapping$schmidtea_gene %in% b_overlap,] 

        target_orth_names=mapping$human_gene
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))  
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task17_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/cat_to_tiger.csv',header=T)
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]  
        a_overlap=intersect(rownames(input_a_rds),mapping$cat_gene)
        input_a_rds=input_a_rds[mapping$cat_gene,]  

        mapping=mapping[mapping$cat_gene %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping$tiger_gene)
        input_b_rds=input_b_rds[mapping$tiger_gene,] 
        mapping=mapping[mapping$tiger_gene %in% b_overlap,] 

        target_orth_names=mapping$cat_gene
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)  
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))  
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task18_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/cat_to_dog.csv',header=T)
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]   
        a_overlap=intersect(rownames(input_a_rds),mapping$cat_gene)
        input_a_rds=input_a_rds[mapping$cat_gene,]  

        mapping=mapping[mapping$cat_gene %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping$dog_gene)
        input_b_rds=input_b_rds[mapping$dog_gene,] 
        mapping=mapping[mapping$dog_gene %in% b_overlap,] 

        target_orth_names=mapping$cat_gene
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task19_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/human_to_MF.csv',header=T)
        mapping=mapping[!duplicated(mapping[,c(4,5)]),] 
        a_name='human_gene'
        b_name='MF_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task20_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/human_to_MM.csv',header=T)
        mapping=mapping[!duplicated(mapping[,c(4,5)]),] 
        a_name='human_gene'
        b_name='MM_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task21_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/human_to_mouse_mm10.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='human_gene'
        b_name='mouse_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 
        mapping=mapping[!duplicated(mapping[,b_name]),]
        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task22_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/mouse_to_MF.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='MF_gene'
        b_name='mouse_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 
        mapping=mapping[!duplicated(mapping[,b_name]),]  

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)  
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task23_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/mouse_to_MM.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='MM_gene'
        b_name='mouse_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 
        mapping=mapping[!duplicated(mapping[,b_name]),]
        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task24_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/human_to_pig.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='human_gene'
        b_name='pig_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task25_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/mouse_to_pig.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='mouse_gene'
        b_name='pig_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task26_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/pig_to_MF.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='MF_gene'
        b_name='pig_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task27_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/pig_to_MM.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='MM_gene'
        b_name='pig_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)  
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task28_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/human_to_fish.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='fish_gene'
        b_name='human_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task29_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/fish_to_MF.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='fish_gene'
        b_name='MF_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    }  else if (target %in% c('task30_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/fish_to_cat.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='cat_gene'
        b_name='fish_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task31_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/pig_to_fish.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='fish_gene'
        b_name='pig_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } else if (target %in% c('task32_')){
        mapping=read.csv('../OrthoFinder/one2one_orthologs/fish_to_MM.csv',header=T)
        #mapping=mapping[!duplicated(mapping),] 
        mapping=mapping[!duplicated(mapping[,c(4,5)]),]
        a_name='fish_gene'
        b_name='MM_gene'  
        a_overlap=intersect(rownames(input_a_rds),mapping[,a_name])
        input_a_rds=input_a_rds[mapping[,a_name],]  

        mapping=mapping[mapping[,a_name] %in% a_overlap,]  
        b_overlap=intersect(rownames(input_b_rds),mapping[,b_name])
        input_b_rds=input_b_rds[mapping[,b_name],] 
        mapping=mapping[mapping[,b_name] %in% b_overlap,] 

        target_orth_names=mapping[,a_name]
        input_b_rds=RenameGenesSeurat(input_b_rds,newnames=target_orth_names)   
        
        common_genes=intersect(rownames(input_a_rds),rownames(input_b_rds))   
        input_a_rds=input_a_rds[common_genes,]
        input_b_rds=input_b_rds[common_genes,]

        input_a_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_a_rds[["RNA"]]))
        input_b_rds@assays$RNA@meta.features <- data.frame(row.names = rownames(input_b_rds[["RNA"]]))
    } 


    ## run the integration
    target_list=list(input_a_rds,input_b_rds)
    integrated_list=seurat4_preprocess(target_list,task_name=gsub("_","",target),batch_label=batch_label,celltype_label=celltype_label,
                                        save_prefix='../output/method_outputs/seurat4/')
    print('integration finished!')
    

}
