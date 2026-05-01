##### the code is for the final visualization for evaluation results of all methods
library(ggplot2)
library(data.table)
library("ggrepel")   

### step 1: combine all the evaluation matrix for all the methods together
setwd('../output/final_evaluation_table/')
target_files=list.files(pattern = '*final_task3_38.csv')
target_files

# Read each file as a data frame and add a label column based on the file name
ldf <- lapply(target_files, function(x) {
  df <- read.csv(x,header=T)
  df$method <- basename(x)
  return(df)
})

# Row bind all the data frames into one
result <- do.call(rbind, ldf)
result$method=gsub('_final_task3_38.csv','',result$method)

dir.create('./processed_final')
write.csv(result,'./processed_final/all_process_result.csv')


### step2: split the total dataframe by task 
dir.create('./visualization')
setwd('./visualization/')
groups <- split(result, result[,'X'])
for (name in names(groups)) {
  write.csv(groups[[name]], paste0(name, '.csv'), row.names = FALSE)
}


### step3: plot the summary matric in csv and png for each task 
# use the code from scib
library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)


source("../../script/scIB_knit_table.R")
# Please put knit_table.R in your working dir
dir.create('./plots')
outdir='./plots'

# get metrics names from columns
metrics <- colnames(metrics_tab_lab)[-1]
metrics <- gsub("\\.", "/", metrics) # replace the . with /
metrics <- gsub("_", " ", metrics)
metrics <- plyr::mapvalues(metrics, from = c("ari batch", "ari celltype", "nmi batch", "nmi celltype", "asw batch", "asw celltype", "lisi batch", "lisi celltype", "kBET", "PCR batch","graph conn","hvg overlap",'trajectory'), 
                           to = c("Batch ARI", "Cell type ARI", "Batch NMI", "Cell type NMI", "Batch ASW", "Cell type ASW", "Graph iLISI", "Graph cLISI","kBET",'PCR batch','Graph connectivity','HVG conservation','Trajectory conservation'))


# metrics names as they are supposed to be ordered
group_batch <- c("Batch ARI", "Batch ASW", "Graph iLISI", "Batch NMI", "kBET","PCR batch","Graph connectivity")
group_bio <- c("Cell type ARI", "Cell type ASW", "Graph cLISI", "Cell type NMI", "HVG conservation",'Trajectory conservation')
# set original values of number of metrics
n_metrics_batch_original <- sum(group_batch %in% metrics)
n_metrics_bio_original <- sum(group_bio %in% metrics)

# order metrics present in the table
matching.order <- match(c(group_batch, group_bio), metrics)
metrics.ord <- metrics[matching.order[!is.na(matching.order)]]

# get methods info from last column
methods_info_full  <- as.character(metrics_tab_lab[,'method'])

# in case methods names start with /
if(substring(methods_info_full[1], 1, 1) == "/"){
  methods_info_full <- sub("/", "", methods_info_full)
}


# data scenarios to be saved in file name
data.scenarios <- unique(metrics_tab_lab$X) # get all the task name
setwd('./plots')
###### Plot one figure for each data task
for (dt.sc in data.scenarios){
  #ind.scen <- grep(paste0("^",dt.sc,'$'), methods_info_full)
  #methods_info <- methods_info_full[ind.scen]
  ind.scen = metrics_tab_lab$X==dt.sc
  metrics_tab_sub <- metrics_tab_lab[ind.scen, ]
  
  methods <- metrics_tab_sub$method
  
  methods_name <- metrics_tab_sub$method
  methods_name <- capitalize(methods_name)
  methods_name <- plyr::mapvalues(methods_name, 
                                  from = c("Bbknn","Fastmnn","Saturn","ScGen","ScVI","Seurat4",'Harmony','Samap'), 
                                  to = c("BBKNN","fastMNN","SATURN",'scGen','scVI',"Seurat v4 CCA","Harmony","SAMap"))
  
  
  ##### Create dataframe 
  metrics_tab <- as.data.frame(metrics_tab_sub[, -1])
  metrics_tab[metrics_tab == ""] <- NA
  colnames(metrics_tab) <- metrics
  metrics_tab=metrics_tab[,-which(colnames(metrics_tab) == "method")]
  #add Methods column
  metrics_tab <- add_column(metrics_tab, "Method" = methods_name, .before = 1)
  
  # reorder columns by metrics
  col.ordered <- c("Method", metrics.ord)
  metrics_tab <- metrics_tab[, col.ordered]
  
  ## Remove columns that are full NAs
  na.col <- apply(metrics_tab, 2, function(x) sum(is.na(x)) == nrow(metrics_tab))
  # redefine numbers of metrics per group
  if(sum(colnames(metrics_tab)[na.col] %in%  group_batch) > 0){
    n_metrics_batch <- n_metrics_batch_original - sum(colnames(metrics_tab)[na.col] %in%  group_batch)
  } else {
    n_metrics_batch <- n_metrics_batch_original
  }
  
  if(sum(colnames(metrics_tab)[na.col] %in% group_bio) > 0){
    n_metrics_bio <- n_metrics_bio_original - sum(colnames(metrics_tab)[na.col] %in% group_bio)
  } else{
    n_metrics_bio <- n_metrics_bio_original
  }
  
  metrics_tab <- metrics_tab[, !na.col]
  
  
  ## Scores should be already scaled [0,1] - however, we aim to compute the scores based on the min-max scaled metrics
  scaled_metrics_tab <- as.matrix(metrics_tab[, -1])
  scaled_metrics_tab <- apply(scaled_metrics_tab, 2, function(x) scale_minmax(x))
  
  # calculate average score by group and overall
  score_group_batch <- rowMeans(scaled_metrics_tab[, 1:n_metrics_batch], na.rm = T)
  score_group_bio <- rowMeans(scaled_metrics_tab[, (1+n_metrics_batch):ncol(scaled_metrics_tab)], 
                              na.rm = T)
  
  weight_batch=0.4
  score_all <- (weight_batch*score_group_batch + (1-weight_batch)*score_group_bio)
  
  metrics_tab <- add_column(metrics_tab, "Overall Score" = score_all, .after = "Method")
  metrics_tab <- add_column(metrics_tab, "Batch Correction" = score_group_batch, .after = "Overall Score")
  metrics_tab <- add_column(metrics_tab, "Bio conservation" = score_group_bio, .after = 3+n_metrics_batch)
  
  # order methods by the overall score
  metrics_tab <- metrics_tab[order(metrics_tab$`Overall Score`,  decreasing = T), ]
  
  write.csv(metrics_tab, file = paste0(outdir, "/", dt.sc, "_summary_scores.csv"), quote = F)
  
  # Delete rows that are empty
  rowsNA <- which(is.na(metrics_tab$`Overall Score`))
  if(length(rowsNA) >0){
    metrics_tab <- metrics_tab[-rowsNA, ]
  }
  
  
  # Defining column_info, row_info and palettes
  row_info <- data.frame(id = metrics_tab$Method)
  
  column_info <- data.frame(id = colnames(metrics_tab),
                            group = c("Text", "Score overall", 
                                      rep("Removal of batch effects", (1 + n_metrics_batch)),
                                      rep("Cell type label variance", (1 + n_metrics_bio))), 
                            geom = c("text", "bar", "bar", 
                                     rep("circle", n_metrics_batch), "bar", rep("circle", n_metrics_bio)),
                            width = c(6,3,3, rep(1,n_metrics_batch), 3, rep(1,n_metrics_bio)),
                            overlay = F)
  
  # defining colors palette
  palettes <- list("Score overall" = "YlGnBu",
                   "Removal of batch effects" = "BuPu",
                   "Cell type label variance" = "RdPu")
  
  row.names(metrics_tab) <- 1:nrow(metrics_tab)
  #source('scIB_knit_table.R')  # need to change soure f
  g <- scIB_knit_table(data = metrics_tab, column_info = column_info, row_info = row_info, palettes = palettes, usability = F)  
  g
  now <- Sys.time()
  #ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dt.sc, "_summary_metrics.pdf"), g, device = cairo_pdf, width = 297, height = 420, units = "mm")
  #ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dt.sc, "_summary_metrics.tiff"), g, device = "tiff", dpi = "retina", width = 297, height = 420, units = "mm")
  ggsave(paste0(outdir,"/", dt.sc, "_summary_metrics.png"), g, device = "png", dpi = "retina", width = 297, height = 420, units = "mm",bg = 'white')
  
  
}



### step 4: plot the figure for each category (cross-family, cross-order, cross-class, cross-genus, cross-phylum)
################################
plot_average_score_table<-function(
  input_table,
  table_save_path,
  figure_save_path){
  colnames(input_table)=c('Method','Overall Score','Batch correction',"Batch ARI", "Batch ASW", "Graph iLISI", "Batch NMI", "kBET",
                          "PCR batch","Graph connectivity",'Bio conservation',"Cell type ARI", "Cell type ASW", "Graph cLISI", "Cell type NMI", "HVG conservation",'Trajectory conservation')
  
  # order methods by the overall score
  input_table <- input_table[order(input_table$`Overall Score`,  decreasing = T), ]
  
  write.csv(input_table, file = table_save_path, quote = F)
  # plot the table plot of the cross-genus result
  # metrics names as they are supposed to be ordered
  group_batch <- c("Batch ARI", "Batch ASW", "Graph iLISI", "Batch NMI", "kBET","PCR batch","Graph connectivity")
  group_bio <- c("Cell type ARI", "Cell type ASW", "Graph cLISI", "Cell type NMI", "HVG conservation",'Trajectory conservation')
  # set original values of number of metrics
  metrics <- colnames(input_table)[-1]
  n_metrics_batch_original <- sum(group_batch %in% metrics)
  n_metrics_bio_original <- sum(group_bio %in% metrics)
  
  methods <- input_table$method
  na.col <- apply(input_table, 2, function(x) sum(is.na(x)) == nrow(input_table))
  # redefine numbers of metrics per group
  if(sum(colnames(input_table)[na.col] %in%  group_batch) > 0){
    n_metrics_batch <- n_metrics_batch_original - sum(colnames(input_table)[na.col] %in%  group_batch)
  } else {
    n_metrics_batch <- n_metrics_batch_original
  }
  
  if(sum(colnames(input_table)[na.col] %in% group_bio) > 0){
    n_metrics_bio <- n_metrics_bio_original - sum(colnames(input_table)[na.col] %in% group_bio)
  } else{
    n_metrics_bio <- n_metrics_bio_original
  }
  metrics_tab=input_table
  metrics_tab <- metrics_tab[, !na.col]
  # Delete rows that are empty
  rowsNA <- which(is.na(metrics_tab$`Overall Score`))
  if(length(rowsNA) >0){
    metrics_tab <- metrics_tab[-rowsNA, ]
  }
  
  # Defining column_info, row_info and palettes
  row_info <- data.frame(id = metrics_tab$Method)
  
  column_info <- data.frame(id = colnames(metrics_tab),
                            group = c("Text", "Score overall", 
                                      rep("Removal of batch effects", (1 + n_metrics_batch)),
                                      rep("Cell type label variance", (1 + n_metrics_bio))), 
                            geom = c("text", "bar", "bar", 
                                     rep("circle", n_metrics_batch), "bar", rep("circle", n_metrics_bio)),
                            width = c(6,3,3, rep(1,n_metrics_batch), 3, rep(1,n_metrics_bio)),
                            overlay = F)
  
  # defining colors palette
  palettes <- list("Score overall" = "YlGnBu",
                   "Removal of batch effects" = "BuPu",
                   "Cell type label variance" = "RdPu")
  
  row.names(metrics_tab) <- 1:nrow(metrics_tab)
  source('../scIB_knit_table.R')
  g <- scIB_knit_table(data = metrics_tab, column_info = column_info, row_info = row_info, palettes = palettes, usability = F)  
  g
  ggsave(figure_save_path, g, device = "png", dpi = "retina", width = 297, height = 420, units = "mm",bg = 'white')
  return(g)
}

# read the data of the cross genus
setwd('../')
dir.create('./plots/summary_plot/',recursive = T)

task3=read.csv('./plots/task3_summary_scores.csv',header = T,row.names = 1)
task17=read.csv('./plots/task17_summary_scores.csv',header = T,row.names = 1)

# merge the data
genus=rbind(task3,task17)
genus_avg=aggregate(genus[, 2:17], list(genus$Method), mean)

plot_average_score_table(genus_avg,table_save_path = './plots/summary_plot/cross_genus_summary_scores.csv',
                         figure_save_path = './plots/summary_plot/cross_genus_avg_summary_metrics.png')
# the figure is the average performance in all the cross-genus tasks, which is Figure 2a

## plot the cross family
task7=read.csv('./plots/task7_summary_scores.csv',header = T,row.names = 1)
task18=read.csv('./plots/task18_summary_scores.csv',header = T,row.names = 1)
task19=read.csv('./plots/task19_summary_scores.csv',header = T,row.names = 1)
task20=read.csv('./task20_summary_scores.csv',header = T,row.names = 1)

family=rbind(task7,task18,task19,task20)
family_avg=aggregate(family[, 2:17], list(family$Method), mean)
plot_average_score_table(family_avg,table_save_path = './plots/summary_plot/cross_family_summary_scores.csv',
                         figure_save_path = './plots/summary_plot/cross_family_avg_summary_metrics.png')
##the figure saved is Figure 3a

summary_family <- family %>%
  group_by(Method) %>%
  dplyr::summarize(
    OverallMean=mean(Overall.Score),
    OverallSD=sd(Overall.Score)/sqrt(n()),
    BatchMean = mean(Batch.Correction),
    BioMean = mean(Bio.conservation),
    BatchSD = sd(Batch.Correction) / sqrt(n()),
    BioSD = sd(Bio.conservation) / sqrt(n())
  )

library('tidyr')
mean_test=summary_family %>% pivot_longer(cols = c('OverallMean','BatchMean','BioMean'),names_to = 'type',values_to = 'mean_value')
sd_test=summary_family %>% pivot_longer(cols = c('OverallSD','BatchSD','BioSD'),names_to = 'type',values_to = 'sd_value')
mean_test=as.data.frame(mean_test[,c(1,5,6)])
sd_test=as.data.frame(sd_test[,c(1,5,6)])
mean_test$type=gsub('OverallMean','Overall score',mean_test$type)
mean_test$type=gsub('BioMean','Bio-conservation score',mean_test$type)
mean_test$type=gsub('BatchMean','Batch correction score',mean_test$type)
sd_test$type=gsub('OverallSD','Overall score',sd_test$type)
sd_test$type=gsub('BioSD','Bio-conservation score',sd_test$type)
sd_test$type=gsub('BatchSD','Batch correction score',sd_test$type)
mean_test$id=paste0(mean_test$Method,'_',mean_test$type)
sd_test$id=paste0(sd_test$Method,'_',sd_test$type)
family_summary=merge(mean_test,sd_test,by='id')

family_summary=family_summary[,c(2,3,4,7)]
head(family_summary)
family_summary$type.x=gsub(' score','',family_summary$type.x)
family_summary$type.x=gsub('Bio-conservation','Bio conservation',family_summary$type.x)

family_overall=mean(family_avg$Overall.Score)
family_batch=mean(family_avg$Batch.Correction)
family_bio=mean(family_avg$Bio.conservation)

family_summary$Method.x=factor(family_summary$Method.x,levels = c('Harmony','SAMap','BBKNN','fastMNN','Scanorama','scVI','Seurat v4 CCA','scGen','SATURN'))
ggplot(family_summary, aes(x = factor(Method.x), y = mean_value, fill = type.x)) + 
  geom_bar(stat = "identity", position = "dodge")  +
  geom_errorbar(aes(ymin=mean_value-sd_value, ymax=mean_value+sd_value), position = position_dodge(0.9), width = 0.25,
                show.legend = FALSE,color='brown') +
  labs(x="Method", y="Score") +
  geom_hline(yintercept = family_batch, linetype = "dashed", color = "#8c96c7",size=1.2,show.legend = T, aes(linetype = "Mean")) +
  geom_hline(yintercept = family_bio, linetype = "dashed", color = "#f868a1",size=1.2,show.legend = T, aes(linetype = "Mean")) +
  geom_hline(yintercept = family_overall, linetype = "dashed", color = "#ffA500",size=1.2,show.legend = T, aes(linetype = "Mean")) +
  scale_fill_manual(labels=c('Batch correction','Bio conservation','Overall'),values = c( "#8c96c7","#f868a1",'#ffA500')) +
  theme_bw() + 
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    legend.position  = "top",
    axis.title       = element_text(size = 28),
    axis.text        = element_text(size = 26),
    axis.title.x = element_blank(),
    panel.border     = element_rect(fill = NA),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(size = 10, colour = "white"),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 18)
  )

ggsave('./plots/summary_plot/cross_family_bar_plot_for_batch_bio_all.png',
       width = 810, height = 297, units = "mm",dpi=300,bg = "white")
### this is the supplementary figure 2e

### plot the cross-order
task6=read.csv('./plots/task6_summary_scores.csv',header = T,row.names = 1)
task21=read.csv('./plots/task21_summary_scores.csv',header = T,row.names = 1)
task22=read.csv('./plots/task22_summary_scores.csv',header = T,row.names = 1)
task23=read.csv('./plots/task23_summary_scores.csv',header = T,row.names = 1)
task24=read.csv('./plots/task24_summary_scores.csv',header = T,row.names = 1)
task25=read.csv('./plots/task25_summary_scores.csv',header = T,row.names = 1)
task26=read.csv('./plots/task26_summary_scores.csv',header = T,row.names = 1)
task27=read.csv('./plots/task27_summary_scores.csv',header = T,row.names = 1)
task4_1=read.csv('./plots/task4_summary_scores.csv',header = T,row.names = 1)

df=rbind(task4_1,task6,task21,task22,task23,task24,task25,task26,task27)
df_avg=aggregate(df[, 2:17], list(df$Method), mean)
plot_average_score_table(df_avg,table_save_path = './plots/summary_plot/cross_order_summary_scores.csv',
                         figure_save_path = './plots/summary_plot/cross_order_avg_summary_metrics.png')
## the figure is Supplementary figure 3b


# plot the scatter plot of the average scores with error bars
# Calculate average and standard error
library(dplyr)
library(plyr)
summary_data <- df %>%
  group_by(Method) %>%
  dplyr::summarize(
    BatchMean = mean(Batch.Correction),
    BioMean = mean(Bio.conservation),
    BatchSD = sd(Batch.Correction) / sqrt(n()),
    BioSD = sd(Bio.conservation) / sqrt(n())
  )
mean_bio=mean(summary_data$BioMean)
mean_batch=mean(summary_data$BatchMean)
# Plotting
library(ggplot2)
library(paletteer)
ggplot(summary_data) +
  aes(
    x      = BatchMean,
    y      = BioMean,
    colour = Method
  ) +
  geom_errorbarh(
    aes(xmin = BatchMean - BatchSD, xmax = BatchMean + BatchSD)
  ) +
  geom_errorbar(
    aes(ymin = BioMean   - BioSD, ymax = BioMean   + BioSD)
  ) +
  geom_point(size = 3, stroke = 1, fill = "white") +
  paletteer::scale_color_paletteer_d(
    palette = "ggsci::category20_d3",
    labels = summary_data$Method
  ) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "Batch correction",
    y = "Bio conservation"
  ) +
  guides(
    colour = guide_legend(
      title          = "Method",
      title.position = "top",
      ncol           = 2,
      order          = 10
    )
  ) +
  theme_minimal() +
  geom_hline(yintercept = mean_bio, linetype = "dashed", color = "red",size=0.6,show.legend = T, aes(linetype = "Mean")) +
  geom_vline(xintercept = mean_batch, linetype = "dashed", color = "red",size=0.6,show.legend = T, aes(linetype = "Mean")) +
  scale_linetype_manual(values = c("dashed"), guide = guide_legend(override.aes = list(color = c("red")))) +
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    legend.position  = "bottom",
    axis.title       = element_text(size = 24),
    axis.text        = element_text(size = 16),
    panel.border     = element_rect(fill = NA),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(size = 16, colour = "white"),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 16)
  )

ggsave('./plots/summary_plot/cross_order_avg_scatter_plot.png',
       width = 210, height = 297, units = "mm",dpi=300,bg = "white")
### the figure saved is figure 3b


##### plot the cross class
#task8=read.csv('./plots/task8_summary_scores.csv',header = T,row.names = 1)
task8=read.csv('./plots/task8_summary_scores.csv',header = T,row.names = 1)
task28=read.csv('./plots/task28_summary_scores.csv',header = T,row.names = 1)
task29=read.csv('./plots/task29_summary_scores.csv',header = T,row.names = 1)
task30=read.csv('./plots/task30_summary_scores.csv',header = T,row.names = 1)
task31=read.csv('./plots/task31_summary_scores.csv',header = T,row.names = 1)
task32=read.csv('./plots/task32_summary_scores.csv',header = T,row.names = 1)
df=rbind(task8,task28,task29,task30,task31,task32)
df_avg=aggregate(df[, 2:17], list(df$Method), mean)
plot_average_score_table(df_avg,table_save_path = './plots/summary_plot/cross_class_summary_scores.csv',
                         figure_save_path = './plots/summary_plot/cross_class_avg_summary_metrics.png')
# the figure saved is supplementary Figure 4b

df_avg=df_avg[order(df_avg$Batch.Correction,decreasing = T),]
df_avg$Group.1=gsub('Seurat v4 CCA','Seurat v4',df_avg$Group.1)
# draw the bos plot of the cross class score
# Plotting boxplots for both scores
library(data.table)
df_long=melt(setDT(df[,c(1,3,11)]),id.vars = c('Method'),variable.name = 'Type')
df_long<-as.data.frame(df_long)
df_long$Method=gsub('Seurat v4 CCA','Seurat v4',df_long$Method)
df_long$Method<-factor(df_long$Method,levels=df_avg$Group.1)
dodge <- position_dodge(width = 0.9)
mean_bio=mean(df_avg$Bio.conservation)
mean_batch=mean(df_avg$Batch.Correction)

ggplot(df_long, aes(x = Method,y=value)) +
  #geom_boxplot(aes(fill = reorder(Type,value)), position = dodge) +
  geom_boxplot(aes(fill = Type), position = dodge) +
  labs( y = "Score") +
  scale_fill_manual(labels=c('Batch correction','Bio conservation'),values = c( "#8c96c7","#f868a1")) +
  theme_minimal() +
  theme(legend.position = "top")+
  guides(
    colour = guide_legend(
      title          = "Method",
      title.position = "top",
      ncol           = 2,
      order          = 10
    )
  ) +
  geom_hline(yintercept = mean_bio, linetype = "dashed", color = "#f868a1",size=1.5,show.legend = T, ) +
  geom_hline(yintercept = mean_batch, linetype = "dashed", color = "#8c96c7",size=1.5,show.legend = T, ) +
  scale_linetype_manual(values = c("dashed",'dashed'),labels=c('avg_batch','avg_bio'),
                        guide = guide_legend(override.aes = list(color = c("#f868a1","#8c96c7")))) +
  theme(
    panel.background = element_rect(colour = 'black',size = 1),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    legend.position  = "top",
    axis.title       = element_text(size = 36),
    axis.text        = element_text(size = 32),
    #axis.text.x        = element_text(face = 'bold'),
    axis.title.x        = element_blank(),
    panel.border     = element_rect(fill = NA),
    strip.background = element_rect(fill = "black"),
    strip.text       = element_text(size = 10, colour = "white"),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 30),
    axis.ticks = element_line(color = 'black',size = 0.5)
  )+
  theme(axis.ticks.length.x =unit(0.15, "cm"),
        #axis.ticks.length.y =unit(-0.25, "cm"),
        legend.key.size = unit(3,"line")
  )
ggsave('./plots/summary_plot/cross_class_avg_box_plot.png',
       width = 579, height = 410, units = "mm",dpi=300,bg = "white")
# the figure saved is figure 3c


###### plot the cross phylum 
task10=read.csv('./plots/task10_summary_scores.csv',header = T,row.names = 1)
task9_1=read.csv('./plots/task9-1_summary_scores.csv',header = T,row.names = 1)
task13=read.csv('./plots/task13_summary_scores.csv',header = T,row.names = 1)
task33=read.csv('./plots/task33_summary_scores.csv',header = T,row.names = 1)
task34=read.csv('./plots/task34_summary_scores.csv',header = T,row.names = 1)
task35=read.csv('./plots/task35_summary_scores.csv',header = T,row.names = 1)
task36=read.csv('./plots/task36_summary_scores.csv',header = T,row.names = 1)
task37=read.csv('./plots/task37_summary_scores.csv',header = T,row.names = 1)
task38=read.csv('./plots/task38_summary_scores.csv',header = T,row.names = 1)
task9=read.csv('./plots/task9_summary_scores.csv',header = T,row.names = 1)
task14=read.csv('./plots/task14_summary_scores.csv',header = T,row.names = 1)
task15=read.csv('./plots/task15_summary_scores.csv',header = T,row.names = 1)
task16=read.csv('./plots/task16_summary_scores.csv',header = T,row.names = 1)
task12=read.csv('./plots/task12_summary_scores.csv',header = T,row.names = 1)
task11=read.csv('./plots/task11_summary_scores.csv',header = T,row.names = 1)

df=rbind(task13,task33,task34,task35,task36,task37,task38,task9,task9_1,task10,task14,task15,task16,task12,task11)
df_avg=aggregate(df[, 2:17], list(df$Method), mean)
plot_average_score_table(df_avg,table_save_path = './plots/summary_plot/cross_phylum_summary_scores.csv',
                         figure_save_path = './plots/summary_plot/cross_phylum_avg_summary_metrics.png')
# the figure saved is supplementary figure 5a

### plot the lollipop plot for all the methods
# Calculate average and standard error
summary_data <- df %>%
  group_by(Method) %>%
  dplyr::summarize(
    BatchMean = mean(Batch.Correction),
    BioMean = mean(Bio.conservation),
    BatchSD = sd(Batch.Correction) / sqrt(n()),
    BioSD = sd(Bio.conservation) / sqrt(n())
  )
mean_bio=mean(summary_data$BioMean)
mean_batch=mean(summary_data$BatchMean)
median_bio=median(summary_data$BioMean)
median_batch=median(summary_data$BatchMean)

library(tidyverse)
library(ggpubr)
axis_margin <- 5.5
df_avg=df_avg[order(df_avg$Batch.Correction,decreasing = F),]
summary_data$Method=factor(summary_data$Method,levels = df_avg$Group.1)

p1<-ggplot(summary_data, aes(BatchMean,Method)) +
  geom_col(fill="#8c96c7",width = 0.5) +
  geom_errorbarh(aes(xmin = BatchMean-BatchSD,xmax = BatchMean+BatchSD), height = 0.2)+
  scale_x_reverse(limits=c(0.8,0)) + #limits=c(0.8,0),expand=c(0,0)
  scale_y_discrete(position = "right") +
  theme(
    #axis.text.y = element_blank(),
    axis.text.y.right  = element_text(margin = margin(0, 2, 0, axis_margin),size=14),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 28),
    plot.margin = margin(axis_margin, 0, axis_margin, axis_margin),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

p2 <- ggplot(summary_data, aes(BioMean,Method)) +
  geom_col(fill='#f868a1',width = 0.5) +
  geom_errorbarh(aes(xmin = BioMean-BioSD,xmax = BioMean+BioSD), height = 0.2)+
  xlim(0,0.8)+
  scale_y_discrete(position = "left") +
  #coord_cartesian(xlim = c(0, 0.8))+
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 28),
    plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
    axis.text.y.left = element_text(margin = margin(0, 2, 0, axis_margin),size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )

ggarrange(p1, p2,ncol = 2)

ggsave('./plots/summary_plot/cross_phyum_bar_plot_batch_bio.png',
       width = 579, height = 410, units = "mm",dpi=300,bg = "white")
# the figure doesn't get shown in the paper since i choose the lollopop plot

## plot the lollipop plot 
p3<-ggplot(summary_data,aes(x = Method, y = BatchMean)) +
  geom_segment(aes(xend = Method,yend=0)) +
  geom_point(size = 6, color = "orange") +
  scale_x_reverse()+
  theme_bw() +
  xlab("") +
  coord_flip() +
  #scale_y_continuous(expand = c(0,0.001))+
  scale_x_discrete(position = "top") +
  scale_y_reverse()+
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 28),
    plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
    axis.text.y.right = element_text(margin = margin(0, 2, 0, axis_margin),size=14),
    panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )+
  theme(panel.border = element_blank(), axis.line = element_line())+
  geom_hline(yintercept = mean_batch,linetype='dashed',color='red',size=1)
#geom_hline(yintercept = median_batch,linetype='dashed',color='green',size=1)

p4<-ggplot(summary_data,aes(x = Method, y = BioMean)) +
  geom_segment(aes(xend = Method,yend=0)) +
  geom_point(size = 6, color = "orange") +
  #scale_y_continuous(expand=c(0,0.001))+
  theme_bw() +
  xlab("") +
  coord_flip() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 28),
    #axis.text.y = element_text(size = 28),
    plot.margin = margin(axis_margin, axis_margin, axis_margin, 0),
    axis.text.y.left = element_text(margin = margin(0, 2, 0, axis_margin),size=14),
    panel.grid.major = element_blank(), #panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")
  )+
  theme(panel.border = element_blank(), axis.line = element_line())+
  geom_hline(yintercept = mean_bio,linetype='dashed',color='red',size=1)
#geom_hline(yintercept = median_bio,linetype='dashed',color='green',size=1)

ggarrange(p3,p4,ncol = 2)
ggsave('./plots/summary_plot/cross_phyum_lopllipop_plot_batch_bio.png',
       width = 400, height = 250, units = "mm",dpi=300,bg = "white")
## the figure saved is figure 4a

### plot the bar plot for the cross-phylum plot (figure 4b)
#total=list(task13,task33,task34,task35,task36,task37,task9,task9_1,task10,task14,task15,task16,task12,task11)
total=list(task19,task21,task24,task28,task14,task15,task16)

factor_levels=c("SATURN","SAMap","fastMNN","Scanorama","BBKNN","scVI","scGen","Seurat v4 CCA","Harmony")

total<-lapply(total, function(df){
  df$Method=factor(df$Method,levels=factor_levels)
  df[, 3] <- rowMeans(df[, c(4:10)], na.rm = TRUE)  
  df[,c(11)]=rowMeans(df[,c(12:17)],na.rm = T)
  df=df[,c(1,3,11)]
  #df$Bio.conservation <- 0.6*df$Bio.conservation
  #df$Batch.Correction <- 0.4*df$Batch.Correction
  return(df)
})

#names(total)=c('task13','task33','task34','task35','task36','task37','task9','task9_1','task10','task14','task15','task16','task12','task11')
names(total)=c('task19','task21','task24','task28','task14','task15','task16')

library(dplyr)
library(tidyr)
# Convert dataframe to long format
#lapply(total, function(df) {
for (i in seq_along(total)) {
  df <- total[[i]]
  df_long <- gather(df, key = score_type, value = score, Batch.Correction, Bio.conservation)
  
  ggplot(df_long) +
    geom_bar(aes(x = Method, y = score, fill = score_type), position = "dodge", stat = "identity") +
    scale_fill_manual(values=c( "#8c96c7","#f868a1"))+
    #scale_fill_manual(values = c("Batch.Correction" = "#8c96c7", "Bio.conservation" = "#f868a1"), name = "Score Type") +
    labs(x = "Method", y = "Overall score") +
    theme_minimal()+
    theme(
      panel.background = element_blank(),  # Remove panel background
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      axis.line = element_line(color = "black"),  # Set axis lines to black color
      axis.text = element_text(size = 19,colour = 'black'),  # Adjust axis text size
      axis.title = element_text(size = 20),  # Adjust axis title size
      axis.ticks = element_line(color = "black"),  # Add tick bars
      axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
      axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
      legend.title = element_text(size = 21),  # Adjust legend title size
      legend.text = element_text(size = 19),  # Adjust legend text size
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.ticks.x  = element_blank(),
    )+
    guides(fill=guide_legend(title=""))+
    scale_fill_discrete(labels=c('Batch correction','Bio conservation'))+
    ylim(c(0,1))+
    labs(title=NULL,x=NULL,y=NULL)+ #title = names(total)[i]
    theme(axis.text.x = element_blank())+
    scale_fill_manual(values=c( "#8c96c7","#f868a1"))
  ggsave(paste0('./plots/summary_plot/cross_phylum_',names(total)[i],'_bar_overall.png'),dpi=300,width = 6299,height = 1024,units = 'px')
  
}


############################################
## below is the dot plots for different tasks
############################################
data=read.csv('./plots/task6_summary_scores.csv',header = T,row.names = 1)
head(data)
data=data[,c(1,3,11)]

# draw the dot plot 
data$Method <- factor(data$Method, levels = unique(data$Method))

library(ggrepel)
# Create the dot plot for task 6: nematostella and stylophora
ggplot(data, aes(x = Batch.Correction, y = Bio.conservation, color = Method)) +
  geom_point(size = 5) +
  labs(x = "Batch correction", y = "Bio-conservation", color = "Method") +
  scale_color_brewer(palette = "Paired") +
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 14),  # Adjust axis text size
    axis.title = element_text(size = 20),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_text(size = 21),  # Adjust legend title size
    legend.text = element_text(size = 19)  # Adjust legend text size
  ) +
  #geom_text(aes(label = Method), vjust = -0.90,show.legend = FALSE,size=5) + # Display method names above the dots
  geom_text_repel(aes(label = Method), size = 7, box.padding = 0.5,show.legend = FALSE) +
  guides(color = guide_legend(override.aes = list(shape = 16)))  # Use pure dots in the legend

ggsave('./plots/summary_plot/task6_overall_dotplot.png',dpi=300)
# the plot doesn't get shown in the paper.

# for task 8: frog and zebrafish
data=read.csv('./plots/task8_summary_scores.csv',header = T,row.names = 1)
head(data)
data
data$Method=factor(data$Method,levels = data$Method)
data=data[,c(1,3,11)]

library(tidyr)
# Convert dataframe to long format
df_long <- gather(data, key = score_type, value = score, Batch.Correction, Bio.conservation)

ggplot(df_long) +
  geom_bar(aes(x = Method, y = score, fill = score_type), position = "dodge", stat = "identity") +
  #scale_fill_manual(values = c("Batch.Correction" = "green", "Bio.conservation" = "blue"), name = "Score Type") +
  labs(x = "Method", y = "Overall score") +
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 14,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 20),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_text(size = 21),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  guides(fill=guide_legend(title=""))+
  scale_fill_discrete(labels=c('Batch correction','Bio-conservation'))

ggsave('./plots/summary_plot/task8_overall_barplot.png',dpi=300)
### the plot doesn't get shown in the paper


# plot the linechart 
task13=read.csv('./plots/task13_summary_scores.csv',header = T,row.names = 1)
task10=read.csv('./plots/task10_summary_scores.csv',header = T,row.names = 1)
task14=read.csv('./plots/task14_summary_scores.csv',header = T,row.names = 1)
task15=read.csv('./plots/task15_summary_scores.csv',header = T,row.names = 1)
task16=read.csv('./plots/task16_summary_scores.csv',header = T,row.names = 1)
task12=read.csv('./plots/task12_summary_scores.csv',header = T,row.names = 1)
task11=read.csv('./plots/task11_summary_scores.csv',header = T,row.names = 1)


task3=read.csv('./plots/task3_summary_scores.csv',header = T,row.names = 1)
task7=read.csv('./plots/task7_summary_scores.csv',header = T,row.names = 1)
task6=read.csv('./plots/task6_summary_scores.csv',header = T,row.names = 1)
task8=read.csv('./plots/task8_summary_scores.csv',header = T,row.names = 1)

task9=read.csv('./plots/task9_summary_scores.csv',header = T,row.names = 1)
task9_1=read.csv('./plots/task9-1_summary_scores.csv',header = T,row.names = 1)

task19=read.csv('./plots/task19_summary_scores.csv',header = T,row.names = 1)
task20=read.csv('./plots/task20_summary_scores.csv',header = T,row.names = 1)
task21=read.csv('./plots/task21_summary_scores.csv',header = T,row.names = 1)
task24=read.csv('./plots/task24_summary_scores.csv',header = T,row.names = 1)
task28=read.csv('./plots/task28_summary_scores.csv',header = T,row.names = 1)

factor_levels=c("SATURN","SAMap","scGen","Seurat v4 CCA","Scanorama","scVI","fastMNN","BBKNN","Harmony")

all_human_inte=list(task19,task21,task24,task28,task14,task15,task16,task11)
total=list(task13,task10,task14,task15,task16,task12,task11)
total_low=list(task3,task7,task5,task6,task8)
total_quality=list(task9,task9_1)

names(all_human_inte)=c('task19','task21','task24','task28','task14','task15','task16','task11')
names(total)=c('task13','task10','task14','task15','task16','task12','task11')
names(total_low)=c('task3','task7','task5','task6','task8')
names(total_quality)=c('task9','task9_1')

for (i in seq_along(all_human_inte)) {
  df <- all_human_inte[[i]]
  df$Method=factor(df$Method,levels=factor_levels)
  df=df[,c(1,11)]
  df$task=names(all_human_inte)[i]
  all_human_inte[[i]]=df
}

library(dplyr)
combined_df <- bind_rows(all_human_inte)
head(combined_df)
colnames(combined_df)[2]<-'Overall.Score'
#combined_df$task<-factor(combined_df$task,levels = c('task13','task10','task14','task15','task16','task12','task11'))
#combined_df$task<-factor(combined_df$task,levels = c('task3','task7','task5','task6','task8'))
#combined_df$task<-factor(combined_df$task,levels = c('task9','task9_1'))
combined_df$task<-factor(combined_df$task,levels = c('task19','task20','task21','task24','task28','task14','task15','task16','task11'))

# Create a line chart for figure 4c and 4d, supplementary 5b
#combined_df=combined_df[combined_df$task %in% c('task14','task15','task16','task11'),]
ggplot(combined_df, aes(x = task, y = Overall.Score, group = Method, color = Method)) +
  geom_line(size=1) +
  geom_point(size=3)+
  labs(x = "Task", y = "Bio-conservation score") +
  scale_color_discrete(name = "Method")+
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 26,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 28),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_blank(),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    axis.text.x = element_text(angle = 45, hjust = 1),
  )+
  guides(fill=guide_legend(title=""))+
  scale_x_discrete(expand = c(0,0.1))+
  theme(legend.position = "top")
ggsave(paste0('./plots/summary_plot/cross_phylum_all_human_line_plot_overall_bioconservation.png'),dpi=300,
       width = 500, height = 200, units = "mm",bg = "white")

total_quality=list(task9,task9_1)
names(total_quality)=c('task9','task9_1')

for (i in seq_along(total_quality)) {
  df <- total_quality[[i]]
  df$Method=factor(df$Method,levels=factor_levels)
  df=df[,c(1,11)]
  df$task=names(total_quality)[i]
  total_quality[[i]]=df
}

combined_df <- bind_rows(total_quality)
head(combined_df)
colnames(combined_df)[2]<-'Overall.Score'
combined_df$task<-factor(combined_df$task,levels = c('task9','task9_1'))
combined_df$Method<-factor(combined_df$Method,levels = c('BBKNN','fastMNN','Harmony','SAMap','SATURN','Scanorama','scGen','scVI','Seurat v4 CCA'))

# Create a line chart for supplementary figure 7 efg
ggplot(combined_df, aes(x = task, y = Overall.Score, group = Method, color = Method)) +
  geom_line(size=1) +
  geom_point(size = 4)+
  scale_color_brewer(palette="Set1")+
  labs(x = "Task", y = "Bio-conservation score") +
  #scale_color_discrete(name = "Method")+
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 26,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 28),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_blank(),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )+
  guides(fill=guide_legend(title=""))+
  scale_x_discrete(expand = c(0,0.1))+
  theme(legend.position = "right")
ggsave(paste0('./plots/summary_plot/data_quality_line_plot_overall_bio.png'),dpi=300,
       width = 3544,height = 2500, units = "px",bg = "white")


total_quality=list(task9,task9_1)
names(total_quality)=c('task9','task9_1')

for (i in seq_along(total_quality)) {
  df <- total_quality[[i]]
  df$Method=factor(df$Method,levels=factor_levels)
  df=df[,c(1,3)]
  df$task=names(total_quality)[i]
  total_quality[[i]]=df
}

combined_df <- bind_rows(total_quality)
head(combined_df)
colnames(combined_df)[2]<-'Overall.Score'
combined_df$task<-factor(combined_df$task,levels = c('task9','task9_1'))
combined_df$Method<-factor(combined_df$Method,levels = c('BBKNN','fastMNN','Harmony','SAMap','SATURN','Scanorama','scGen','scVI','Seurat v4 CCA'))

# Create a line chart for supplementary figure 7 efg
ggplot(combined_df, aes(x = task, y = Overall.Score, group = Method, color = Method)) +
  geom_line(size=1) +
  geom_point(size = 4)+
  scale_color_brewer(palette="Set1")+
  labs(x = "Task", y = "Batch correction score") +
  #scale_color_discrete(name = "Method")+
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 26,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 28),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_blank(),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )+
  guides(fill=guide_legend(title=""))+
  scale_x_discrete(expand = c(0,0.1))+
  theme(legend.position = "right")
ggsave(paste0('./plots/summary_plot/data_quality_line_plot_overall_batch.png'),dpi=300,
       width = 3544,height = 2500, units = "px",bg = "white")


total_quality=list(task9,task9_1)
names(total_quality)=c('task9','task9_1')

for (i in seq_along(total_quality)) {
  df <- total_quality[[i]]
  df$Method=factor(df$Method,levels=factor_levels)
  df=df[,c(1,2)]
  df$task=names(total_quality)[i]
  total_quality[[i]]=df
}

combined_df <- bind_rows(total_quality)
head(combined_df)
colnames(combined_df)[2]<-'Overall.Score'
combined_df$task<-factor(combined_df$task,levels = c('task9','task9_1'))
combined_df$Method<-factor(combined_df$Method,levels = c('BBKNN','fastMNN','Harmony','SAMap','SATURN','Scanorama','scGen','scVI','Seurat v4 CCA'))

# Create a line chart for supplementary figure 7 efg
ggplot(combined_df, aes(x = task, y = Overall.Score, group = Method, color = Method)) +
  geom_line(size=1) +
  geom_point(size = 4)+
  scale_color_brewer(palette="Set1")+
  labs(x = "Task", y = "Overall score") +
  #scale_color_discrete(name = "Method")+
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 26,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 28),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_blank(),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )+
  guides(fill=guide_legend(title=""))+
  scale_x_discrete(expand = c(0,0.1))+
  theme(legend.position = "right")
ggsave(paste0('./plots/summary_plot/data_quality_line_plot_overall.png'),dpi=300,
       width = 3544,height = 2500, units = "px",bg = "white")

# fro task4
task4=read.csv('./plots/task4_summary_scores.csv',header = T,row.names = 1)
head(task4)

data=task4
data$Method=factor(data$Method,levels = data$Method)
data=data[,c(1,3,11)]

#data$Bio.conservation <- 0.6*data$Bio.conservation
#data$Batch.Correction <- 0.4*data$Batch.Correction

library(tidyr)
# Convert dataframe to long format
df_long <- gather(data, key = score_type, value = score, Batch.Correction, Bio.conservation)
df_long$Method
ggplot(df_long) +
  geom_bar(aes(x = Method, y = score, fill = score_type), position = "dodge", stat = "identity") +
  scale_fill_manual(values = c("Batch.Correction" = "#8c96c7","Bio.conservation" = "#f868a1"), name = "Score Type") +
  labs(x = "Method", y = "Overall score") +
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 26,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 28),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_text(size = 21),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )+
  guides(fill=guide_legend(title=""))
#scale_fill_discrete(labels=c('Batch correction','Bio-conservation'))

ggsave('./plots/summary_plot/task4_overall.png',dpi=300,
       width = 4512,height = 2800, units = "px")
# the figure saved is Figure 5b

## draw the unbalanced changes line plot of all methods, figure 5d
task13=read.csv('./plots/task13_summary_scores.csv',header = T,row.names = 1)
task33=read.csv('./plots/task33_summary_scores.csv',header = T,row.names = 1)
task34=read.csv('./plots/task34_summary_scores.csv',header = T,row.names = 1)
task35=read.csv('./plots/task35_summary_scores.csv',header = T,row.names = 1)
task36=read.csv('./plots/task36_summary_scores.csv',header = T,row.names = 1)
task37=read.csv('./plots/task37_summary_scores.csv',header = T,row.names = 1)
task38=read.csv('./plots/task38_summary_scores.csv',header = T,row.names = 1)

task13=task13[,c(1,2,3,11)]
task33=task33[,c(1,2,3,11)]
task34=task34[,c(1,2,3,11)]
task35=task35[,c(1,2,3,11)]
task36=task36[,c(1,2,3,11)]
task37=task37[,c(1,2,3,11)]
task38=task38[,c(1,2,3,11)]

# combined them together
library(dplyr)
# totally fish has 1143079 cells, urchin has 60399 cells
# the experiment ratio is to two species cell ratio
merged_unbalamced_data <- bind_rows(
  mutate(task33, Experiment = 1), #1%
  mutate(task34, Experiment = 6), #6%
  mutate(task35, Experiment = 20), #20%
  mutate(task36, Experiment = 40), #40%
  mutate(task37, Experiment = 60), #60%
  mutate(task38, Experiment = 80) #80%
)

merged_unbalamced_data$Method=factor(merged_unbalamced_data$Method)
# plot the overall score of the data 
ggplot(merged_unbalamced_data, aes(x = Experiment, y = Bio.conservation,group=Method,color=Method)) +
  geom_line(size=1) +
  geom_point(size = 4)+
  scale_color_brewer(palette="Set1")+
  #scale_x_continuous(breaks = c(0.16, 0.99, 1.1, 3.3, 6.6, 9.9, 13.2)) +
  scale_x_continuous(breaks = c(1,6,20,40,60, 80)) +
  geom_point() +  # Add points for each data point
  labs(x = "Subsample percentage(%)", y = "Bio-conservation score") +
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 26,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 28),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_blank(),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    #axis.text.x = element_text(angle = 0, hjust = 1),
  )+
  guides(fill=guide_legend(title="")) +
  theme(legend.position = "top")
ggsave('./plots/summary_plot/unbalance_dataset_bio-conservation_overall.png',
       dpi=300,width = 3544,height = 2500,units = 'px')
# the saved figure is figure 5d

ggplot(merged_unbalamced_data, aes(x = Experiment, y = Batch.Correction,group=Method,color=Method)) +
  geom_line(size=1) +
  geom_point(size = 4)+
  scale_color_brewer(palette="Set1")+
  #scale_x_continuous(breaks = c(0.16, 0.99, 1.1, 3.3, 6.6, 9.9, 13.2)) +
  scale_x_continuous(breaks = c(1,6,20,40,60, 80)) +
  geom_point() +  # Add points for each data point
  labs(x = "Subsample percentage(%)", y = "Batch correction score") +
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 26,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 28),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_blank(),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    #axis.text.x = element_text(angle = 0, hjust = 1),
  )+
  guides(fill=guide_legend(title="")) +
  theme(legend.position = "top")
ggsave('./plots/summary_plot/unbalance_dataset_batch_correction_overall.png',
       dpi=300,width = 3544,height = 2500,units = 'px')
# the saved figure is supplementary figure 7d

ggplot(merged_unbalamced_data, aes(x = Experiment, y = Overall.Score,group=Method,color=Method)) +
  geom_line(size=1) +
  geom_point(size = 4)+
  scale_color_brewer(palette="Set1")+
  #scale_x_continuous(breaks = c(0.16, 0.99, 1.1, 3.3, 6.6, 9.9, 13.2)) +
  scale_x_continuous(breaks = c(1,6,20,40,60, 80)) +
  geom_point() +  # Add points for each data point
  labs(x = "Subsample percentage(%)", y = "Overall score") +
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 26,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 28),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_blank(),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    #axis.text.x = element_text(angle = 0, hjust = 1),
  )+
  guides(fill=guide_legend(title="")) +
  theme(legend.position = "top")
ggsave('./plots/summary_plot/unbalance_dataset_overall.png',
       dpi=300,width = 3544,height = 2500,units = 'px')
# the saved figure is supplementary figure 7c

# draw the overall score barplot for all the tasks
task13=read.csv('./plots/task13_summary_scores.csv',header = T,row.names = 1)
task33=read.csv('./plots/task33_summary_scores.csv',header = T,row.names = 1)
task34=read.csv('./plots/task34_summary_scores.csv',header = T,row.names = 1)
task35=read.csv('./plots/task35_summary_scores.csv',header = T,row.names = 1)
task36=read.csv('./plots/task36_summary_scores.csv',header = T,row.names = 1)
task37=read.csv('./plots/task37_summary_scores.csv',header = T,row.names = 1)
task38=read.csv('./plots/task38_summary_scores.csv',header = T,row.names = 1)
task9=read.csv('./plots/task9_summary_scores.csv',header = T,row.names = 1)
task14=read.csv('./plots/task14_summary_scores.csv',header = T,row.names = 1)
task15=read.csv('./plots/task15_summary_scores.csv',header = T,row.names = 1)
task16=read.csv('./plots/task16_summary_scores.csv',header = T,row.names = 1)
task12=read.csv('./plots/task12_summary_scores.csv',header = T,row.names = 1)
task11=read.csv('./plots/task11_summary_scores.csv',header = T,row.names = 1)

task3=read.csv('./plots/task3_summary_scores.csv',header = T,row.names = 1)
task17=read.csv('./plots/task17_summary_scores.csv',header = T,row.names = 1)

task7=read.csv('./plots/task7_summary_scores.csv',header = T,row.names = 1)
task18=read.csv('./plots/task18_summary_scores.csv',header = T,row.names = 1)
task19=read.csv('./plots/task19_summary_scores.csv',header = T,row.names = 1)
task20=read.csv('./plots/task20_summary_scores.csv',header = T,row.names = 1)

task6=read.csv('./plots/task6_summary_scores.csv',header = T,row.names = 1)
task21=read.csv('./plots/task21_summary_scores.csv',header = T,row.names = 1)
task22=read.csv('./plots/task22_summary_scores.csv',header = T,row.names = 1)
task23=read.csv('./plots/task23_summary_scores.csv',header = T,row.names = 1)
task24=read.csv('./plots/task24_summary_scores.csv',header = T,row.names = 1)
task25=read.csv('./plots/task25_summary_scores.csv',header = T,row.names = 1)
task26=read.csv('./plots/task26_summary_scores.csv',header = T,row.names = 1)
task27=read.csv('./plots/task27_summary_scores.csv',header = T,row.names = 1)

task8=read.csv('./plots/task8_summary_scores.csv',header = T,row.names = 1)
task10=read.csv('./plots/task10_summary_scores.csv',header = T,row.names = 1)
task4_1=read.csv('./plots/task4_summary_scores.csv',header = T,row.names = 1)
task9_1=read.csv('./plots/task9-1_summary_scores.csv',header = T,row.names = 1)
task28=read.csv('./plots/task28_summary_scores.csv',header = T,row.names = 1)
task29=read.csv('./plots/task29_summary_scores.csv',header = T,row.names = 1)
task30=read.csv('./plots/task30_summary_scores.csv',header = T,row.names = 1)
task31=read.csv('./plots/task31_summary_scores.csv',header = T,row.names = 1)
task32=read.csv('./plots/task32_summary_scores.csv',header = T,row.names = 1)


average_function=function(
  other_task_list,
  atlas_task_list,
  target_column_name
){
  other_tasks <- do.call(rbind, other_task_list)
  atlas_tasks <- do.call(rbind,atlas_task_list)
  other_tasks2 <- aggregate(get(target_column_name) ~ Method, data = other_tasks, FUN = mean)
  atlas_tasks2 <- aggregate(get(target_column_name) ~ Method, data = atlas_tasks, FUN = mean)
  all_tasks <- list(other_tasks2,atlas_tasks2)
  all_tasks <- do.call(rbind, all_tasks)
  colnames(all_tasks)[2]=target_column_name
  final <- aggregate(get(target_column_name) ~ Method, data = all_tasks, FUN = mean)
  colnames(final)[2]=target_column_name
  return(final)
}

genus_list <- list(task3,task17)  # Replace with your actual dataframes
genus <- do.call(rbind, genus_list)
genus <- aggregate(Overall.Score ~ Method, data = genus, FUN = mean)

family_list <- list(task7,task18,task19,task20)  # Replace with your actual dataframes
family <- do.call(rbind, family_list)
family <- aggregate(Overall.Score ~ Method, data = family, FUN = mean)

family=average_function(other_task_list = list(task18,task19,task20),
                      atlas_task_list = list(task7),
                      target_column_name = 'Overall.Score')

order_list <- list(task4_1,task6,task21,task22,task23,task24,task25,task26,task27)  # Replace with your actual dataframes
order <- do.call(rbind, order_list)
order <- aggregate(Overall.Score ~ Method, data = order, FUN = mean)

order=average_function(other_task_list = list(task21,task22,task23,task24,task25,task26,task27),
                        atlas_task_list = list(task4_1,task6),
                        target_column_name = 'Overall.Score')

class_list <- list(task8,task28,task29,task30,task31,task32)  # Replace with your actual dataframes
class <- do.call(rbind, class_list)
class <- aggregate(Overall.Score ~ Method, data = class, FUN = mean)

class=average_function(other_task_list = list(task28,task29,task30,task31,task32),
                       atlas_task_list = list(task8),
                       target_column_name = 'Overall.Score')

phylum_list <- list(task13,task33,task34,task35,task36,task37,task38,task9,task9_1,task10,task14,task15,task16,task12,task11)  # Replace with your actual dataframes
phylum <- do.call(rbind, phylum_list)
phylum <- aggregate(Overall.Score ~ Method, data = phylum, FUN = mean)

phylum=average_function(other_task_list = list(task15,task9,task9_1),
                       atlas_task_list = list(task13,task33,task34,task35,task36,task37,task38,task10,task14,task15,task16,task12,task11),
                       target_column_name = 'Overall.Score')

#dataframes <- list(task3, task7, task6, task8,task13,task9,task9_1,task10,task14,task4_1,task15,task16,task12,task11)  # Replace with your actual dataframes
dataframes <- list(genus,family,order,class,phylum)

# Row bind all the data frames into one
result <- do.call(rbind, dataframes)

averages <- aggregate(Overall.Score ~ Method, data = result, FUN = mean)
errors <- aggregate(Overall.Score ~ Method, data = result, FUN = function(x) sd(x) / sqrt(length(x)))
results <- merge(averages, errors, by = "Method")

library(ggplot2)

results$Method<-factor(results$Method,levels = results[order(results$Overall.Score.x,decreasing = T),][,1])

ggplot(results, aes(x = Method, y = Overall.Score.x, color = Method)) +
  geom_bar(stat = "identity", position = "dodge",fill='white',size=1,width = 0.8) +
  geom_errorbar(aes(ymin = Overall.Score.x - Overall.Score.y, ymax = Overall.Score.x + Overall.Score.y), width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Method", y = "Overall score", fill = "Method") +
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 14,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 20),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    legend.title = element_text(size = 21),  # Adjust legend title size
    legend.text = element_text(size = 19),  # Adjust legend text size
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank()
  )+
  guides(fill=guide_legend(title=""))+
  theme(legend.position = "none")
ggsave('./plots/summary_plot/barplot_with_error_bar_overall_all_methods.png',dpi=300,width = 3562, height = 2870, units = "px",bg = 'white')
# the figure doesn't get shown in the paper, since it is the overll scores 

# plot the average bioconservation and batch correction
genus <- do.call(rbind, genus_list)
genus <- aggregate(Batch.Correction ~ Method, data = genus, FUN = mean)

family <- do.call(rbind, family_list)
family <- aggregate(Batch.Correction ~ Method, data = family, FUN = mean)

order <- do.call(rbind, order_list)
order <- aggregate(Batch.Correction ~ Method, data = order, FUN = mean)

class <- do.call(rbind, class_list)
class <- aggregate(Batch.Correction ~ Method, data = class, FUN = mean)

phylum <- do.call(rbind, phylum_list)
phylum <- aggregate(Batch.Correction ~ Method, data = phylum, FUN = mean)

family=average_function(other_task_list = list(task18,task19,task20),
                        atlas_task_list = list(task7),
                        target_column_name = 'Batch.Correction')

order=average_function(other_task_list = list(task21,task22,task23,task24,task25,task26,task27),
                       atlas_task_list = list(task4_1,task6),
                       target_column_name = 'Batch.Correction')

class=average_function(other_task_list = list(task28,task29,task30,task31,task32),
                       atlas_task_list = list(task8),
                       target_column_name = 'Batch.Correction')

phylum=average_function(other_task_list = list(task15,task9,task9_1),
                        atlas_task_list = list(task13,task33,task34,task35,task36,task37,task38,task10,task14,task15,task16,task12,task11),
                        target_column_name = 'Batch.Correction')

dataframes <- list(genus,family,order,class,phylum)

result <- do.call(rbind, dataframes)

errors_batch <- aggregate(Batch.Correction ~ Method, data = result, FUN = function(x) sd(x) / sqrt(length(x)))
averages_batch <- aggregate(Batch.Correction ~ Method, data = result, FUN = mean)
median_batch=median(result$Batch.Correction)
mean_batch=mean(result$Batch.Correction)

genus <- do.call(rbind, genus_list)
genus <- aggregate(Bio.conservation ~ Method, data = genus, FUN = mean)

family <- do.call(rbind, family_list)
family <- aggregate(Bio.conservation ~ Method, data = family, FUN = mean)

order <- do.call(rbind, order_list)
order <- aggregate(Bio.conservation ~ Method, data = order, FUN = mean)

class <- do.call(rbind, class_list)
class <- aggregate(Bio.conservation ~ Method, data = class, FUN = mean)

phylum <- do.call(rbind, phylum_list)
phylum <- aggregate(Bio.conservation ~ Method, data = phylum, FUN = mean)

family=average_function(other_task_list = list(task18,task19,task20),
                        atlas_task_list = list(task7),
                        target_column_name = 'Bio.conservation')

order=average_function(other_task_list = list(task21,task22,task23,task24,task25,task26,task27),
                       atlas_task_list = list(task4_1,task6),
                       target_column_name = 'Bio.conservation')

class=average_function(other_task_list = list(task28,task29,task30,task31,task32),
                       atlas_task_list = list(task8),
                       target_column_name = 'Bio.conservation')

phylum=average_function(other_task_list = list(task15,task9,task9_1),
                        atlas_task_list = list(task13,task33,task34,task35,task36,task37,task38,task10,task14,task15,task16,task12,task11),
                        target_column_name = 'Bio.conservation')

dataframes <- list(genus,family,order,class,phylum)

result <- do.call(rbind, dataframes)
errors_bio <- aggregate(Bio.conservation ~ Method, data = result, FUN = function(x) sd(x) / sqrt(length(x)))
averages_bio <- aggregate(Bio.conservation ~ Method, data = result, FUN = mean)
median_bio=median(result$Bio.conservation)
mean_bio=mean(result$Bio.conservation)

results <- merge(averages_batch, averages_bio, by = "Method")
results

results$Method=factor(results$Method,levels = c('SATURN','SAMap','Harmony','Scanorama','BBKNN','scGen','Seurat v4 CCA','scVI','fastMNN'))

ggplot(results, aes(x = Batch.Correction, y = Bio.conservation, color = Method)) +
  geom_point(size = 5) +
  #geom_errorbar(aes(ymin = Bio.conservation - errors_bio, ymax = Bio.conservation + errors_bio), width = 0) +
  #geom_errorbarh(aes(xmin = Batch.Correction - errors_batch, xmax = Batch.Correction + errors_batch), height=0) +
  labs(x = "Batch correction score", y = "Bio-conservation score", color = "Method") +
  geom_hline(yintercept = mean_bio, linetype = "dashed", color = "red",size=1,show.legend = T, aes(linetype = "Mean")) +
  geom_vline(xintercept = mean_batch, linetype = "dashed", color = "red",size=1,show.legend = T, aes(linetype = "Mean")) +
  scale_linetype_manual(values = c("dashed"), guide = guide_legend(override.aes = list(color = c("red")))) +
  scale_color_brewer(palette = "Paired") +
  theme_minimal()+
  theme(
    panel.background = element_blank(),  # Remove panel background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    axis.line = element_line(color = "black"),  # Set axis lines to black color
    axis.text = element_text(size = 24,colour = 'black'),  # Adjust axis text size
    axis.title = element_text(size = 28,colour = 'black'),  # Adjust axis title size
    axis.ticks = element_line(color = "black"),  # Add tick bars
    axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
    axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
    #legend.title = element_text(size = 21),  # Adjust legend title size
    #legend.text = element_text(size = 19),  # Adjust legend text size
    panel.border = element_rect(colour = "black", fill=NA, size=2) # add the frame of the plot
  ) +
  theme(legend.position="none")+
  ylim(0,1)+
  xlim(0,1)+
  #geom_text(aes(label = Method), vjust = -0.90,show.legend = FALSE,size=5) + # Display method names above the dots
  geom_text_repel(aes(label = Method), size = 11, box.padding = 0.5,show.legend = FALSE) +
  guides(color = guide_legend(override.aes = list(shape = 16)))  # Use pure dots in the legend

ggsave('./plots/summary_plot/Overall_mean_score_for_all_methods_inall_tasks_atlas_different.png',dpi=300,width = 3245,height = 3045,units = 'px')
# the figure saved is figure 7a

library(dplyr)
dataframes <- list(genus_list,family_list,order_list,class_list,phylum_list)
names(dataframes)=c('cross_genus','cross_family','cross_order','cross_class','cross_phylum')
for (i in 1:length(dataframes)){
  result=Reduce(full_join,dataframes[[i]])
  mean_bio=mean(result$Bio.conservation)
  mean_batch=mean(result$Batch.Correction)
  print(names(dataframes)[i])
  print(mean_bio)
  print(mean_batch)
  result_new=result %>%
    group_by(Method) %>%
    summarise_at(vars(Batch.Correction,Bio.conservation), list(name = mean))
  result_new=as.data.frame(result_new)
  colnames(result_new)[2:3]=c('Batch.Correction','Bio.conservation')
  result_new$Method=factor(result_new$Method,levels = c('SATURN','SAMap','Harmony','Scanorama','BBKNN','scGen','Seurat v4 CCA','scVI','fastMNN'))
  ggplot(result_new, aes(x = Batch.Correction, y = Bio.conservation, color = Method)) +
    geom_point(size = 5) +
    labs(x = "Batch correction score", y = "Bio-conservation score", color = "Method") +
    geom_hline(yintercept = mean_bio, linetype = "dashed", color = "red",size=1,show.legend = T, aes(linetype = "Median")) +
    geom_vline(xintercept = mean_batch, linetype = "dashed", color = "red",size=1,show.legend = T, aes(linetype = "Median")) +
    scale_color_brewer(palette = "Paired") +
    theme_minimal()+
    theme(
      panel.background = element_blank(),  # Remove panel background
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      axis.line = element_line(color = "black"),  # Set axis lines to black color
      axis.text = element_text(size = 14),  # Adjust axis text size
      axis.title = element_text(size = 20),  # Adjust axis title size
      axis.ticks = element_line(color = "black"),  # Add tick bars
      axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
      axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
      legend.title = element_text(size = 21),  # Adjust legend title size
      legend.text = element_text(size = 19),  # Adjust legend text size
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    ) +
    xlim(0,1)+
    ylim(0,1)+
    #geom_text(aes(label = Method), vjust = -0.90,show.legend = FALSE,size=5) + # Display method names above the dots
    geom_text_repel(aes(label = Method), size = 7, box.padding = 0.5,show.legend = FALSE) +
    guides(color = guide_legend(override.aes = list(shape = 16)))  # Use pure dots in the legend
  ggsave(paste0('./plots/summary_plot/supplementary_',names(dataframes)[i],'_overall.png'),dpi=300,width = 3604,height = 2399,units = 'px')
}

dataframes <- list(task3, task17,task7,task18,task19,task20, task6, task21,task22,task23,task24,task25,task26,
                   task27,task8,task10,task4_1,task9_1,task28,task29,task30,task31,task32,task13,task33,task34,
                   task35,task36,task37,task38,task9,task14,task15,task16,task12,task11)  # Replace with your actual dataframes
names(dataframes)=c('task3','task17','task7','task18','task19','task20', 'task6','task21','task22','task23','task24','task25','task26',
                    'task27','task8','task10','task4_1','task9_1','task28','task29','task30','task31','task32','task13','task33','task34',
                    'task35','task36','task37','task38','task9','task14','task15','task16','task12','task11')
for (i in 1:length(dataframes)){
  result=dataframes[[i]]
  mean_bio=mean(result$Bio.conservation)
  mean_batch=mean(result$Batch.Correction)
  print(names(dataframes)[i])
  print(mean_bio)
  print(mean_batch)
  
  result$Method=factor(result$Method,levels = c('SATURN','SAMap','Harmony','Scanorama','BBKNN','scGen','Seurat v4 CCA','scVI','fastMNN'))
  ggplot(result, aes(x = Batch.Correction, y = Bio.conservation, color = Method)) +
    geom_point(size = 5) +
    labs(x = "Batch correction score", y = "Bio-conservation score", color = "Method") +
    geom_hline(yintercept = mean_bio, linetype = "dashed", color = "red",size=1,show.legend = T, aes(linetype = "Median")) +
    geom_vline(xintercept = mean_batch, linetype = "dashed", color = "red",size=1,show.legend = T, aes(linetype = "Median")) +
    scale_color_brewer(palette = "Paired") +
    theme_minimal()+
    theme(
      panel.background = element_blank(),  # Remove panel background
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      axis.line = element_line(color = "black"),  # Set axis lines to black color
      axis.text = element_text(size = 14),  # Adjust axis text size
      axis.title = element_text(size = 20),  # Adjust axis title size
      axis.ticks = element_line(color = "black"),  # Add tick bars
      axis.ticks.length = unit(0.2, "cm"),  # Set tick bar length
      axis.ticks.margin = unit(0.2, "cm"),  # Set tick bar margin
      legend.title = element_text(size = 21),  # Adjust legend title size
      legend.text = element_text(size = 19),  # Adjust legend text size
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    ) +
    xlim(0,1)+
    ylim(0,1)+
    #geom_text(aes(label = Method), vjust = -0.90,show.legend = FALSE,size=5) + # Display method names above the dots
    geom_text_repel(aes(label = Method), size = 7, box.padding = 0.5,show.legend = FALSE) +
    guides(color = guide_legend(override.aes = list(shape = 16)))  # Use pure dots in the legend
  ggsave(paste0('./plots/summary_plot/supplementary_',names(dataframes)[i],'_overall.png'),dpi=300,width = 3604,height = 2399,units = 'px')
}

### plot the overall result table (median overall scores for the genus, family, order, class, phylum)
genus <- do.call(rbind, genus_list)
genus <- aggregate(Overall.Score ~ Method, data = genus, FUN = mean)

family <- do.call(rbind, family_list)
family <- aggregate(Overall.Score ~ Method, data = family, FUN = mean)

order <- do.call(rbind, order_list)
order <- aggregate(Overall.Score ~ Method, data = order, FUN = mean)

class <- do.call(rbind, class_list)
class <- aggregate(Overall.Score ~ Method, data = class, FUN = mean)

phylum <- do.call(rbind, phylum_list)
phylum <- aggregate(Overall.Score ~ Method, data = phylum, FUN = mean)

dataframes <- list(genus,family,order,class,phylum)


colnames(genus)[2]='Cross genus'
colnames(family)[2]='Cross family'
colnames(order)[2]='Cross order'
colnames(class)[2]='Cross class'
colnames(phylum)[2]='Cross phylum'

# Merge dataframes
merged_df <- merge(genus, family, by = "Method", all = TRUE)
merged_df <- merge(merged_df, order, by = "Method", all = TRUE)
merged_df <- merge(merged_df, class, by = "Method", all = TRUE)
result <- merge(merged_df, phylum, by = "Method", all = TRUE)



## process the dataframe
##### Create dataframe 
metrics_tab <- result
metrics_tab$row_sum=rowSums(metrics_tab[,c(2,3,4,5,6)])

# order methods by the overall score
metrics_tab <- metrics_tab[order(metrics_tab$row_sum,  decreasing = T), ]

write.csv(metrics_tab, file = "./plots/summary_plot/all_animal_group_summary_mean_overall_scores.csv", quote = F)
metrics_tab
# Delete row sums
metrics_tab <- metrics_tab[,-c(7) ]

# Defining column_info, row_info and palettes
row_info <- data.frame(id = metrics_tab$Method)

column_info <- data.frame(id = colnames(metrics_tab),
                          group = c("Text", "Cross genus", 
                                    "Cross family", "Cross order", "Cross class", "Cross phylum"), 
                          geom = c("text", "bar", "bar", 
                                   "bar", "bar","bar"),
                          width = c(6,3,3,3,3,3),
                          overlay = F)

# defining colors palette
palettes <- list("Cross genus" = "YlGnBu",
                 "Cross family" = "BuPu",
                 "Cross order" = "RdPu",
                 "Cross class" = "Oranges", 
                 "Cross phylum" = "Greens")

row.names(metrics_tab) <- 1:nrow(metrics_tab)


#source('scIB_knit_table.R')  # need to change soure f
data = metrics_tab
library(dplyr)
library(scales)
library(ggimage)
library(cowplot)


# no point in making these into parameters
row_height <- 1.1
row_space <- .1
row_bigspace <- .5
col_width <- 1.1
col_space <- .2
col_bigspace <- .5
segment_data <- NULL

# DETERMINE ROW POSITIONS
if (!"group" %in% colnames(row_info) || all(is.na(row_info$group))) {
  row_info$group <- ""
  row_groups <- tibble(group = "")
  plot_row_annotation <- FALSE
} else {
  plot_row_annotation <- TRUE
}

row_pos <-   #### need to check further
  row_info %>%
  group_by(group) %>%
  dplyr::mutate(group_i = row_number()) %>% 
  ungroup() %>%
  dplyr::mutate(
    row_i = row_number(),  #
    colour_background = group_i %% 2 == 1,
    do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
    ysep = ifelse(do_spacing, row_height + 2 * row_space, row_space),
    y = - (row_i * row_height + cumsum(ysep)),
    ymin = y - row_height / 2,
    ymax = y + row_height / 2
  )

# DETERMINE COLUMN POSITIONS
if (!"group" %in% colnames(column_info) || all(is.na(column_info$group))) {
  column_info$group <- ""
  plot_column_annotation <- FALSE
} else {
  plot_column_annotation <- TRUE
}

column_info <-
  column_info %>%
  add_column_if_missing(width = col_width, overlay = FALSE)


column_pos <-
  column_info %>%
  mutate(
    do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
    xsep = case_when(
      overlay ~ c(0, -head(width, -1)),
      do_spacing ~ col_bigspace,
      TRUE ~ col_space
    ),
    xwidth = case_when(
      overlay & width < 0 ~ width - xsep,
      overlay ~ -xsep,
      TRUE ~ width
    ),
    xmax = cumsum(xwidth + xsep),
    xmin = xmax - xwidth,
    x = xmin + xwidth / 2
  )


##########################
#### CREATE GEOM DATA ####
##########################

# gather circle data
ind_circle <- which(column_info$geom == "circle")
if(length(ind_circle) > 0){
  dat_mat <- as.matrix(data[, ind_circle])
  col_palette <- data.frame(metric = colnames(dat_mat), 
                            group = column_info[match(colnames(dat_mat), column_info$id), "group"])
  
  col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
  
  circle_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                                  function(x) rep(x, nrow(dat_mat)))), 
                            x0 = unlist(lapply(column_pos$x[ind_circle], 
                                               function(x) rep(x, nrow(dat_mat)))), 
                            y0 = rep(row_pos$y, ncol(dat_mat)),
                            r = row_height/2*as.vector(sqrt(dat_mat))
  )
  for(l in unique(circle_data$label)){
    ind_l <- which(circle_data$label == l)
    circle_data[ind_l, "r"] <- rescale(circle_data[ind_l, "r"], to = c(0.05, 0.55), from = range(circle_data[ind_l, "r"], na.rm = T))
  }
  
  colors <- NULL
  
  
  for(i in 1:ncol(dat_mat)){
    palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
    colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
  }
  
  circle_data$colors <- colors
}


# gather bar data
ind_bar <- which(column_info$geom == "bar")
dat_mat <- as.matrix(data[, ind_bar])

col_palette <- data.frame(metric = colnames(dat_mat), 
                          group = column_info[match(colnames(dat_mat), column_info$id), "group"])

col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])


rect_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                              function(x) rep(x, nrow(dat_mat)))),
                        method = rep(row_info$id, ncol(dat_mat)),
                        value = as.vector(dat_mat),
                        xmin = unlist(lapply(column_pos[ind_bar, "xmin"], 
                                             function(x) rep(x, nrow(dat_mat)))),
                        xmax = unlist(lapply(column_pos[ind_bar, "xmax"], 
                                             function(x) rep(x, nrow(dat_mat)))),
                        ymin = rep(row_pos$ymin, ncol(dat_mat)),
                        ymax = rep(row_pos$ymax, ncol(dat_mat)),
                        xwidth = unlist(lapply(column_pos[ind_bar, "xwidth"], 
                                               function(x) rep(x, nrow(dat_mat))))
)
rect_data <- rect_data %>%
  add_column_if_missing(hjust = 0) %>%
  mutate(
    xmin = xmin + (1 - value) * xwidth * hjust,
    xmax = xmax - (1 - value) * xwidth * (1 - hjust)
  )

colors <- NULL
for(i in 1:ncol(dat_mat)){
  palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
  colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
}

rect_data$colors <- colors



# gather text data
ind_text <- which(column_info$geom == "text")
dat_mat <- as.matrix(data[, ind_text])
atac=FALSE
if(atac){
  colnames(dat_mat)[1] <- "Method"
}
colnames(dat_mat)[1] <- "Method"
text_data <- data.frame(label_value = as.vector(dat_mat), 
                        group = rep(colnames(dat_mat), each = nrow(dat_mat)),
                        xmin = unlist(lapply(column_pos[ind_text, "xmin"], 
                                             function(x) rep(x, nrow(dat_mat)))),
                        xmax = unlist(lapply(column_pos[ind_text, "xmax"], 
                                             function(x) rep(x, nrow(dat_mat)))),
                        ymin = rep(row_pos$ymin, ncol(dat_mat)),
                        ymax = rep(row_pos$ymax, ncol(dat_mat)),
                        size = 5, fontface = "plain", stringsAsFactors = F)

text_data$colors <- "black"

# ADD top3 ranking for each bar column
usability=FALSE
atac_best=FALSE
if(usability || atac_best){
  cols_bar <- unique(rect_data$label)
  cols_bar <- as.character(cols_bar[!is.na(cols_bar)])
  for(c in cols_bar){
    rect_tmp <- rect_data[rect_data$label == c,]
    rect_tmp <- add_column(rect_tmp, "label_value" = as.character(rank(-rect_tmp$value, ties.method = "min")))
    rect_tmp <- rect_tmp[rect_tmp$label_value %in% c("1", "2", "3"), c("label_value", "xmin", "xmax", "ymin", "ymax")]
    rect_tmp <- add_column(rect_tmp, "size" = 2.5, .after = "ymax")
    rect_tmp <- add_column(rect_tmp, "colors" = "black", .after = "size")
    rect_tmp <- add_column(rect_tmp, "fontface" = "plain", .after = "colors")
    rect_tmp <- add_column(rect_tmp, "group" = "top3", .after = "fontface")
    text_data <- bind_rows(text_data, rect_tmp)
  }
}



# ADD COLUMN NAMES
df <- column_pos %>% filter(id != "Method") %>% filter(id != "Ranking")

if (nrow(df) > 0) {
  segment_data <- segment_data %>% bind_rows(
    df %>% transmute(x = x, xend = x, y = -.3, yend = -.1, size = .5)
  )
  text_data <-
    bind_rows(
      text_data,
      df %>% transmute(
        xmin = x, xmax = x, ymin = 0, ymax = -0.5,
        angle = 30, vjust = 0, hjust = 0,
        label_value = id, 
        size = 3
      )
    )
}


# GENERATE ROW ANNOTATION
if (plot_row_annotation) {
  row_annotation <-
    row_pos %>% 
    select(group, ymin, ymax) %>%
    group_by(group) %>%
    summarise(
      ymin = min(ymin),
      ymax = max(ymax),
      y = (ymin + ymax) / 2
    ) %>%
    ungroup() %>%
    mutate(xmin = -.5, xmax = 5) %>%
    filter(!is.na(group), group != "")
  
  text_data <- text_data %>% bind_rows(
    row_annotation %>%
      transmute(xmin, xmax, ymin = ymax + row_space, label_value = group %>% gsub("\n", " ", .), 
                hjust = 0, vjust = .5, fontface = "bold", size = 4) %>%
      mutate(ymax = ymin + row_height)
  )
}

# gather image data
ind_img <- which(column_info$geom == "image")
if(length(ind_img) > 0){
  dat_mat <- as.matrix(data[, ind_img])
  
  image_data <- data.frame(x = unlist(lapply(column_pos$x[ind_img], 
                                             function(x) rep(x, nrow(dat_mat)))), 
                           y = rep(row_pos$y, ncol(dat_mat)),
                           image = mapvalues(dat_mat, from = c("graph", "embed", "gene"), 
                                             to = c("./img/graph.png", "./img/embedding.png", "./img/matrix.png")),
                           stringsAsFactors = FALSE
  )
  
}

suppressWarnings({
  minimum_x <- min(column_pos$xmin, segment_data$x, segment_data$xend, 
                   text_data$xmin, na.rm = TRUE)
  maximum_x <- max(column_pos$xmax, segment_data$x, segment_data$xend, 
                   text_data$xmax, na.rm = TRUE)
  minimum_y <- min(row_pos$ymin, segment_data$y, segment_data$yend,  
                   text_data$ymin, na.rm = TRUE)
  maximum_y <- max(row_pos$ymax, segment_data$y, segment_data$yend, 
                   text_data$ymax, na.rm = TRUE)
})

####################################
###   CREATE HARDCODED LEGENDS   ###
####################################

x_min_output <- minimum_x+0.5
x_min_scaling <- minimum_x + 5.5
x_min_ranking <- ifelse(atac, minimum_x + 5.5, minimum_x + 10.5)
x_min_score <-  ifelse(atac, minimum_x + 11, minimum_x + 17)

leg_max_y <- minimum_y - .5

# Create legend for Output
leg_min_x <- x_min_output
output_title_data <- data.frame(xmin = leg_min_x, 
                                xmax = leg_min_x+ 2, 
                                ymin = leg_max_y - 1, 
                                ymax = leg_max_y, 
                                label_value = "Output", 
                                hjust = 0, vjust = 0, 
                                fontface = "bold",
                                size = 3)

output_img <- data.frame(x = leg_min_x+0.5,
                         y = c(leg_max_y-2, leg_max_y-3.2,leg_max_y-4.4),
                         image = c("./img/matrix.png", "./img/embedding.png", "./img/graph.png")
)
if(atac || atac_best){
  output_text <- data.frame(xmin = leg_min_x+1.5, 
                            xmax = leg_min_x+3, 
                            ymin = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
                            ymax = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
                            label_value = c("feature", "embed", "graph"), 
                            hjust = 0, vjust = 0, 
                            fontface = "plain",
                            size = 3)
} else{
  output_text <- data.frame(xmin = leg_min_x+1.5, 
                            xmax = leg_min_x+3, 
                            ymin = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
                            ymax = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6), 
                            label_value = c("gene", "embed", "graph"), 
                            hjust = 0, vjust = 0, 
                            fontface = "plain",
                            size = 3)
}

text_data <- bind_rows(text_data, output_title_data) #output_text,
#image_data <- bind_rows(image_data, output_img)

# Create legend for scaling
if(!atac && !atac_best){
  leg_min_x <- x_min_scaling
  scaling_title_data <- data.frame(xmin = leg_min_x, 
                                   xmax = leg_min_x+ 2, 
                                   ymin = leg_max_y - 1, 
                                   ymax = leg_max_y, 
                                   label_value = "Scaling", 
                                   hjust = 0, vjust = 0, 
                                   fontface = "bold",
                                   size = 3)
  
  scaling_text <- data.frame(xmin = c(leg_min_x, leg_min_x+1), 
                             xmax = c(leg_min_x+0.5, leg_min_x+3), 
                             ymin = c(rep(leg_max_y-2,2), rep(leg_max_y-3,2)), 
                             ymax = c(rep(leg_max_y-1,2), rep(leg_max_y-2,2)), 
                             label_value = c("+", ": scaled", "-", ": unscaled"), 
                             hjust = 0, vjust = 0, 
                             fontface = c("bold","plain", "bold", "plain"),
                             size = c(5,3,5,3))
  
  #text_data <- bind_rows(text_data, scaling_title_data)  #scaling_text
}

# CREATE LEGEND for ranking colors
leg_min_x <- x_min_ranking
rank_groups <- as.character(column_info[column_info$geom == "bar", "group"])

if(usability){
  rank_minimum_x <- list("RNA" = leg_min_x, 
                         "Simulation" = leg_min_x+1, 
                         "Usability" = leg_min_x+2,
                         "Scalability" = leg_min_x+3)
  leg_max_x <- leg_min_x+3
} else if(atac_best){
  rank_minimum_x <- list("ATAC_windows" = leg_min_x, 
                         "ATAC_peaks" = leg_min_x+1, 
                         "ATAC_genes" = leg_min_x+2)
  leg_max_x <- leg_min_x+2
} else{
  rank_minimum_x <- list("Cross genus" = leg_min_x, 
                         "Cross family" = leg_min_x+1, 
                         "Cross order" = leg_min_x+2,
                         "Cross class" = leg_min_x+3,
                         "Cross phylum" = leg_min_x+4
  )
  leg_max_x <- leg_min_x+4
}

rank_title_data <- data.frame(xmin = leg_min_x, 
                              xmax = leg_min_x+ 2, 
                              ymin = leg_max_y - 1, 
                              ymax = leg_max_y, 
                              label_value = "Ranking", 
                              hjust = 0, vjust = 0, 
                              fontface = "bold")

for(rg in rank_groups){
  rank_palette <- colorRampPalette(rev(brewer.pal(9, palettes[[rg]])))(5)
  
  
  rank_data <- data.frame(xmin = rank_minimum_x[[rg]],
                          xmax = rank_minimum_x[[rg]] + .8,
                          ymin = seq(leg_max_y-4, leg_max_y - 2, by = .5),
                          ymax = seq(leg_max_y-3.5, leg_max_y -1.5, by = .5),
                          border = TRUE,
                          colors = rank_palette
  )
  rect_data <- bind_rows(rect_data, rank_data)
  
}

# create arrow for ranking
arrow_data <- data.frame(x = leg_max_x + 1.5, 
                         xend = leg_max_x +1.5, 
                         y = leg_max_y-4, 
                         yend = leg_max_y -1.5)


# add text next to the arrow
arrow_text <- data.frame(xmin = leg_max_x +2, 
                         xmax = leg_max_x +2.5, 
                         ymin = c(leg_max_y-2, leg_max_y-4), 
                         ymax = c(leg_max_y-1.5, leg_max_y-3.5 ), 
                         label_value = c("1", as.character(nrow(data))), 
                         hjust = 0, vjust = 0, size = 2.5)


text_data <- bind_rows(text_data, rank_title_data, arrow_text)

# CREATE LEGEND for circle scores
# circle legend
if(!usability && !atac_best){
  cir_minimum_x <- x_min_score
  
  cir_legend_size <- 1
  cir_legend_space <- .1
  
  cir_legend_dat <-
    data.frame(
      value = seq(0, 1, by = .2),
      r = row_height/2*seq(0, 1, by = .2)
    )
  cir_legend_dat$r <- rescale(cir_legend_dat$r, to = c(0.05, 0.55), from = range(cir_legend_dat$r, na.rm = T))
  
  x0 <- vector("integer", nrow(cir_legend_dat))
  for(i in 1:length(x0)){
    if(i == 1){
      x0[i] <- cir_minimum_x + cir_legend_space + cir_legend_dat$r[i]
    }
    else {
      x0[i] <- x0[i-1] + cir_legend_dat$r[i-1] + cir_legend_space + cir_legend_dat$r[i]
    }
  }
  
  cir_legend_dat$x0 <- x0
  cir_legend_min_y <- leg_max_y-4
  cir_legend_dat$y0 <- cir_legend_min_y + 1 + cir_legend_dat$r
  
  cir_legend_dat$colors <- NULL
  cir_maximum_x <- max(cir_legend_dat$x0)
  
  cir_title_data <- data_frame(xmin = cir_minimum_x, 
                               xmax = cir_maximum_x, 
                               ymin = leg_max_y -1, 
                               ymax = leg_max_y,
                               label_value = "Score", 
                               hjust = 0, vjust = 0, fontface = "bold")
  
  cir_value_data <- data.frame(xmin = cir_legend_dat$x0 - cir_legend_dat$r,
                               xmax = cir_legend_dat$x0 + cir_legend_dat$r,
                               ymin = cir_legend_min_y,
                               ymax = cir_legend_min_y +3,
                               hjust = .5, vjust = 0, size = 2.5,
                               label_value = ifelse(cir_legend_dat$value %in% c(0, 1), 
                                                    paste0(cir_legend_dat$value*100, "%"), ""))
  
  circle_data <- bind_rows(circle_data, cir_legend_dat)
  text_data <- bind_rows(text_data, cir_title_data, cir_value_data)
  
  
  
}

minimum_y <- min(minimum_y, min(text_data$ymin, na.rm = TRUE))

########################
##### COMPOSE PLOT #####
########################

g <-
  ggplot() +
  coord_equal(expand = FALSE) +
  scale_alpha_identity() +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_size_identity() +
  scale_linetype_identity() +
  cowplot::theme_nothing()

# PLOT ROW BACKGROUNDS
df <- row_pos %>% filter(colour_background)
if (nrow(df) > 0) {
  g <- g + geom_rect(aes(xmin = min(column_pos$xmin)-.25, xmax = max(column_pos$xmax)+.25, ymin = ymin - (row_space / 2), ymax = ymax + (row_space / 2)), df, fill = "#DDDDDD")
} 



# PLOT CIRCLES
if (length(ind_circle) > 0) {
  g <- g + ggforce::geom_circle(aes(x0 = x0, y0 = y0, fill= colors, r = r), circle_data, size=.25)
}


# PLOT RECTANGLES
if (nrow(rect_data) > 0) {
  # add defaults for optional values
  rect_data <- rect_data %>%
    add_column_if_missing(alpha = 1, border = TRUE, border_colour = "black") %>%
    mutate(border_colour = ifelse(border, border_colour, NA))
  
  g <- g + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colors, colour = border_colour, alpha = alpha), rect_data, size = .25)
}


# PLOT TEXT
if (nrow(text_data) > 0) {
  # add defaults for optional values
  text_data <- text_data %>%
    add_column_if_missing(
      hjust = .5,
      vjust = .5,
      size = 3,
      fontface = "plain",
      colors = "black",
      lineheight = 1,
      angle = 0
    ) %>%
    mutate(
      angle2 = angle / 360 * 2 * pi,
      cosa = cos(angle2) %>% round(2),
      sina = sin(angle2) %>% round(2),
      alphax = ifelse(cosa < 0, 1 - hjust, hjust) * abs(cosa) + ifelse(sina > 0, 1 - vjust, vjust) * abs(sina),
      alphay = ifelse(sina < 0, 1 - hjust, hjust) * abs(sina) + ifelse(cosa < 0, 1 - vjust, vjust) * abs(cosa),
      x = (1 - alphax) * xmin + alphax * xmax,
      y = (1 - alphay) * ymin + alphay * ymax
    ) %>%
    filter(label_value != "")
  # Set fontface for legend bold
  text_data[text_data$label_value == "Ranking", "fontface"] <- "bold"
  # Set fontface for ranking numbers bold
  if(usability || atac_best){
    text_data[1:nrow(data), "fontface"] <- "bold"
  }
  # subset text_data to left-aligned rows
  text_data_left <- text_data[which(text_data$group == "Method" | text_data$group == "top3"), ]
  text_data <- text_data[-which(text_data$group == "Method" | text_data$group == "top3"), ]
  
  g <- g + geom_text(aes(x = x, y = y, label = label_value, colour = colors, hjust = hjust, vjust = vjust, size = size, fontface = fontface, angle = angle), data = text_data)
  
  text_data_left[text_data_left$group == "Method", "x"] <- text_data_left[text_data_left$group == "Method", "x"] - 3
  if(usability || atac_best){
    text_data_left[text_data_left$group == "top3", "x"] <- text_data_left[text_data_left$group == "top3", "xmin"] + .3
    text_data_left[text_data_left$group == "Method", "x"] <- text_data_left[text_data_left$group == "Method", "x"] + .5
  }
  g <- g + geom_text(aes(x = x, y = y, label = label_value, colour = colors, hjust = "left", vjust = vjust, size = size, fontface = fontface, angle = angle), data = text_data_left)
}



# PLOT SEGMENTS
if (nrow(segment_data) > 0) {
  # add defaults for optional values
  segment_data <- segment_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")
  
  g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), segment_data)
}

# PLOT ARROW RANKING
if (nrow(arrow_data) > 0) {
  # add defaults for optional values
  arrow_data <- arrow_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")
  
  g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), arrow_data, arrow = arrow(length = unit(0.1, "cm")), lineend = "round", linejoin = "bevel")
}

# PLOT IMAGES
if(length(ind_img) > 0){
  for(r in 1:nrow(image_data)){
    g <- g + cowplot::draw_image(image = image_data$image[r], x = image_data[r, "x"]-.5, y = image_data[r, "y"]-.5)
  }
  
}

# ADD SIZE
# reserve a bit more room for text that wants to go outside the frame
minimum_x <- minimum_x - 2
maximum_x <- maximum_x + 5
minimum_y <- minimum_y - 2
maximum_y <- maximum_y + 4

g$width <- maximum_x - minimum_x
g$height <- maximum_y - minimum_y

g <- g + expand_limits(x = c(minimum_x, maximum_x), y = c(minimum_y, maximum_y))


g

ggsave("./plots/summary_plot/All_summary_mean_overall_score_ranking_table_atlas_different.png", g, device = "png", dpi = "retina", width = 297, height = 420, units = "mm",bg = 'white')
# the figure saved is the figure 7b


#### for radar plot 
genus_list <- list(task3,task17)  # Replace with your actual dataframes
genus <- do.call(rbind, genus_list)
genus$task='cross genus'

family_list <- list(task7,task18,task19,task20)  # Replace with your actual dataframes
family <- do.call(rbind, family_list)
family$task='cross family'

order_list <- list(task4_1,task6,task21,task22,task23,task24,task25,task26,task27)  # Replace with your actual dataframes
order <- do.call(rbind, order_list)
order$task='cross order'

class_list <- list(task8,task28,task29,task30,task31,task32)  # Replace with your actual dataframes
class <- do.call(rbind, class_list)
class$task='cross class'

phylum_list <- list(task13,task33,task34,task35,task36,task37,task38,task9,task9_1,task10,task14,task15,task16,task12,task11)  # Replace with your actual dataframes
phylum <- do.call(rbind, phylum_list)
phylum$task='cross phylum'

#dataframes <- list(task3, task7, task6, task8,task13,task9,task9_1,task10,task14,task4_1,task15,task16,task12,task11)  # Replace with your actual dataframes
dataframes <- list(genus,family,order,class,phylum)
result <- do.call(rbind, dataframes)
result=result[,c(1,2,3,11,18)]
head(result)


bio_result=result[,c(1,4,5)]
bio_result_wide=reshape(bio_result, idvar = "Method", timevar = "task", direction = "wide")
rownames(bio_result_wide)=bio_result_wide$Method
colnames(bio_result_wide)=gsub('Bio.conservation.','',colnames(bio_result_wide))

library(fmsb)
bio_result_wide_norm<- data.frame(sapply(bio_result_wide[,2:6], scales::rescale))
bio_result_wide_norm<-bio_result_wide
rownames(bio_result_wide_norm)=bio_result_wide$Method
bio_result_wide_norm=bio_result_wide_norm[,c(2,3,4,5,6)]
bio_result_wide_norm[nrow(bio_result_wide_norm)+1,]=c(1,1,1,1,1)
bio_result_wide_norm[nrow(bio_result_wide_norm)+1,]=c(0,0,0,0,0)
bio_result_wide_norm=bio_result_wide_norm[c(10,11,1,2,3,4,5,6,7,8,9),]
set1_colors <- brewer.pal(9, "Set1")[1:9]
png('./plots/summary_plot/radar_plot_bio_scores.png',width = 800, height = 800, units = "px")
radarchart(bio_result_wide_norm,
           axistype = 4,
           pcol = set1_colors,
           plwd = 2,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "black",
           frame = FALSE,
           vlcex = 0.8)
legend("topright", legend = rownames(bio_result_wide_norm)[3:11], fill = set1_colors, pch = 16, cex = 1.2, bty = "n")
dev.off()

batch_result=result[,c(1,3,5)]
batch_result_wide=reshape(batch_result, idvar = "Method", timevar = "task", direction = "wide")
rownames(batch_result_wide)=batch_result_wide$Method
colnames(batch_result_wide)=gsub('Batch.Conservation.','',colnames(batch_result_wide))

library(fmsb)
batch_result_wide_norm<- data.frame(sapply(batch_result_wide[,2:6], scales::rescale))
batch_result_wide_norm<-batch_result_wide
rownames(batch_result_wide_norm)=batch_result_wide$Method
batch_result_wide_norm=batch_result_wide_norm[,c(2,3,4,5,6)]
batch_result_wide_norm[nrow(batch_result_wide_norm)+1,]=c(1,1,1,1,1)
batch_result_wide_norm[nrow(batch_result_wide_norm)+1,]=c(0,0,0,0,0)
batch_result_wide_norm=batch_result_wide_norm[c(10,11,1,2,3,4,5,6,7,8,9),]
png('./plots/summary_plot/radar_plot_batch_scores.png',width = 800, height = 800, units = "px")
radarchart(batch_result_wide_norm,
           axistype = 4,
           pcol = set1_colors,
           plwd = 2,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "black",
           frame = FALSE,
           vlcex = 0.8)
legend("topright", legend = rownames(batch_result_wide_norm)[3:11], fill = set1_colors, pch = 16, cex = 1.2, bty = "n")
dev.off()

for (method in unique(result$Method)) {
  print(method)
  result_saturn=result[result$Method==method,]
  result_saturn$task=factor(result_saturn$task)
  result_saturn_ready=result_saturn %>%
    group_by(task) %>%
    dplyr::summarize(mean_batch=mean(Batch.Correction),mean_bio=mean(Bio.conservation))
  result_saturn_t=t(as.data.frame(result_saturn_ready[,-1]))
  colnames(result_saturn_t)=result_saturn_ready$task
  result_saturn_t=as.data.frame(result_saturn_t)
  result_saturn_t[nrow(result_saturn_t)+1,]=c(1,1,1,1,1)
  result_saturn_t[nrow(result_saturn_t)+1,]=c(0,0,0,0,0)
  result_saturn_t=result_saturn_t[c(3,4,'mean_batch','mean_bio'),]
  radarchart(result_saturn_t,
             axistype = 4,
             pcol = 1:4,
             plwd = 2,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "black",
             frame = FALSE,
             vlcex = 0.8)
}




