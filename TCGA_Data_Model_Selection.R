# TCGA Data Model Selection
# A violin diagram fo Class vs Class using the difffernt models and data
# sample1 sample2 data model class1 class2 class_vs_class

setwd("./GliomGEM")
library(devtools)
require(gridExtra)
library(grid)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(knitr)
library(devtools)
library(tidyverse)
library(tidyr)
library(ggplot2)
require(dplyr)
library(rrcov)
library(data.table)
library(vroom)
library(tibble)
library(lemon)
library(ggplotify)
library(patchwork)
#files = list.files('models/',pattern = '_Model_Rxns')
files = list.files('models/',pattern = 'jaccard_matrix')
files <- str_subset(files,'FPKM')
files <- str_subset(files,'Subtypes',negate = TRUE)
# Select 
#patterns <- c('DrugMetabolismRemoved','NoMedium','CSF_Thiele2020'
patterns <- c('NoCuration','CSF_Thiele2020')

Pattern = paste(patterns, collapse="|")
files <- files[grepl(Pattern,files)]


meta = read.csv('data/TCGA_TCGBiolinks_metadata_Summary.csv')
meta$barcode = gsub('\\-','_',meta$barcode)
meta$sample = gsub('\\-','_',meta$sample)
meta$sample = substr(meta$sample ,1,nchar(meta$sample)-1)
meta[meta$TYPE==0,"TYPE"] <- "NA_"
meta[is.na(meta$paper_IDH.status),"paper_IDH.status"] <- "NA_"
meta[meta$paper_X1p.19q.codeletion=="","paper_X1p.19q.codeletion"] <- "NA_"
meta[meta$primary_diagnosis=="","primary_diagnosis"] <- "NA_"
meta[meta$paper_Histology=="","paper_Histology"] <- "NA_"
meta$IDH.status <- meta$paper_IDH.status
meta$X1p.19q.codeletion <- meta$paper_X1p.19q.codeletion

# Create an empty dataframe to append to it
df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(df) <- c("rownames", "colnames",'similarity', "data",'model',"curation")
df <- df %>%
  mutate_if(is.logical, as.character)
df$similarity <- as.numeric(df$similarity )

# Iterate for each model-data reactions and calculate the similarity
for(i in 1:length(files)) {
  # models_rxns <-read.csv(str_c('models/',files[i]),sep=',')
  # rowname= models_rxns[,length(models_rxns)]
  # rownames(models_rxns) = substr(rowname,1,nchar(rowname)-4)
  # models_rxns = models_rxns[,1:length(models_rxns)-1]
  # models_rxns = models_rxns[!rowSums(models_rxns)==0,]
  # #models_rxns = models_rxns[!rowMeans(models_rxns)==1,]
  # jaccard <- function(x, y) { return(length(intersect(x, y)) / length(union(x, y))) }
  # jaccard_sim <- 1-as.matrix(dist(models_rxns, method = "binary") )
  # jaccard_sim <- as.data.frame(jaccard_sim)
  jaccard_sim <-read_csv(str_c('./Sample_models/',files[i]))
  rownames(jaccard_sim) <- jaccard_sim$Row 
  jaccard_sim$rownames <- jaccard_sim$Row 
  jaccard_sim <- jaccard_sim[ , -which(names(jaccard_sim) %in% c("Row"))]
  jaccard_sim %>%
    pivot_longer(!rownames, names_to = "colnames", values_to = "similarity") ->
    jaccard_sim_long
  jaccard_sim_long$similarity <- 1-jaccard_sim_long$similarity
  jaccard_sim_long$rownames <- str_replace(jaccard_sim_long$rownames,'.mat','')
  jaccard_sim_long$colnames <- str_replace(jaccard_sim_long$colnames,'.mat','')
  jaccard_sim_long$data <- str_split(files[i],'_')[[1]][4]
  jaccard_sim_long$model <- str_split(files[i],'_')[[1]][1]
  jaccard_sim_long$curation <- str_replace(str_split(files[i],'_')[[1]][5],"jaccard","")
  df2 <- jaccard_sim_long#[1:floor(nrow(jaccard_sim_long)/2),]
  #df2 = df2[!duplicated(df2$similarity),]
  #df2$rc <- paste(df2$rownames,df2$colnames, sep='')
  #df2$cr <- paste(df2$colnames,df2$rownames, sep='')
  #df2[df2$rc %in% df2$cr,] -> x
  #df2$rc == df2$cr -> y
  #x <- c(df2$rownames == df2$colnames & df2$colnames==df2$rownames)
  #sum(x,na.rm=TRUE)
  df <- bind_rows(df,df2[,1:6])
}

##
unique(df$curation)
max(na.omit(df$similarity))
min(na.omit(df$similarity))
min(df$similarity)
sum(is.na(df$similarity))
df_missing  <- df[is.na(df$similarity),]

## Create heatmaps as subplots
mat_breaks <- seq(0.3,1,length.out = 10)#min(df$similarity), max(df$similarity), length.out = 10)
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

annotation_colors = list( #c('#BB173A','#2F6790','forestgreen')
  WHO_2021 = c(AST_IDH_mut="#2F6790",GBM_IDH_wt='#BB173A',ODG_IDH_mut_Codel="forestgreen",CTRL="#505050",NA_="white"),#,ODG_AST="black"),
  #primary_diagnosis = c(Oligodendroglioma_NOS="aquamarine",Oligodendroglioma_anaplastic="aquamarine4",Astrocytoma_anaplastic = "orange",Astrocytoma_NOS="darkorange",Glioblastoma="red",Solid_Tissue_Normal="green",Mixed_glioma='blue',NA_="white"),
  #paper_Histology = c(astrocytoma = "yellow", glioblastoma="red",oligodendroglioma="green",oligoastrocytoma="black",NA_="white"),
  #TYPE = c(AST="yellow", CTRL="blue",GBM="red",ODG="green",NA_="white",UNK="white"),
  IDH.status = c(Mutant="yellow", WT="red",NA_="white"),
  X1p.19q.codeletion =  c(non_codel="cyan", codel="coral",NA_="white")
)

length(unique(df$rownames))
length(unique(df$colnames))
unique(meta$paper_Histology)
unique(meta$WHO_2021)
save_pheatmap_png <- function(x, filename, units="in", width=10, height=7, res=300) {
  png(filename, width = width, height = height, res = res, units=units)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf <- function(x, filename, units="in", width=10, height=7, res=300) {
  pdf(filename, width = width, height = height, units=units)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
dir.create(file.path('Figure/', 'QC'))

plot_list=list()
files_<- files
## Heatmaps for all glioma classes excluding mixed-glioma
for(i in 1:length(files_)) { #
  # models_rxns <-read.csv(str_c('models/',files[i]),sep=',')
  # rowname= models_rxns[,length(models_rxns)]
  # rownames(models_rxns) = substr(rowname,1,nchar(rowname)-4)
  # models_rxns = models_rxns[,1:length(models_rxns)-1]
  # models_rxns = models_rxns[!rowSums(models_rxns)==0,]
  # #models_rxns = models_rxns[!rowMeans(models_rxns)==1,]
  # jaccard <- function(x, y) { return(length(intersect(x, y)) / length(union(x, y))) }
  # jaccard_sim <- 1-as.matrix(dist(models_rxns, method = "binary") )
  # jaccard_sim <- as.data.frame(jaccard_sim)
  # #jaccard_sim$rownames <- rownames(models_rxns) 
  # #colnames(jaccard_sim) <- substr(colnames(jaccard_sim) ,1,15)
  # #rownames(jaccard_sim) <- substr(rownames(jaccard_sim) ,1,15)
  jaccard_sim <-read_csv(str_c('./Sample_models/',files_[i]))
  files_[i] <- str_replace(files_[i],'NoCuration','NoMedium')
  files_[i] <- str_replace(files_[i],'Biolinks','TCGABiolinks')
  files_[i] <- str_replace(files_[i],'RSEM','Ceccarelli2016')
  
  jaccard_sim <- as.data.frame(jaccard_sim)
  rownames(jaccard_sim) <- jaccard_sim$Row 
  rownames(jaccard_sim) <- str_replace(rownames(jaccard_sim),'.mat','')
  colnames(jaccard_sim) <- str_replace(colnames(jaccard_sim),'.mat','')
  jaccard_sim <- jaccard_sim[ , -which(names(jaccard_sim) %in% c("Row"))]
  jaccard_sim <- 1-jaccard_sim
  #
  ## Selected intersected samples between the jaccard similarity and metadata
    if (nchar(colnames(jaccard_sim)[1])==28) {
      column <- "barcode"
    } else {
      column <- "sample"
    }
  meta_who2021 <- meta[meta$WHO_2021!="NA_",]
  meta_who2021 <- meta_who2021[meta_who2021$WHO_2021!="ODG_AST",] # remove Mixed glioma
  intersected_samples <- intersect(rownames(jaccard_sim),meta_who2021[,column])
  meta_samples <- meta_who2021[meta_who2021[,column] %in% intersected_samples,] 
  
  jaccard_sim_ <- jaccard_sim[rownames(jaccard_sim) %in% intersected_samples, colnames(jaccard_sim) %in% intersected_samples]
  meta_samples<- meta_samples[order(rownames(jaccard_sim_)),]
  rownames(meta_samples) <- rownames(jaccard_sim_)
  meta_samples <- meta_samples[,c('IDH.status','X1p.19q.codeletion','WHO_2021')]#,'OS','Age',,'Study'
  #mat_breaks <- quantile_breaks(jaccard_sim, n = 11)
  #jaccard_sim$rownames <- NULL
  data <- str_split(files_[i],'_')[[1]][4]
  model <- str_split(files_[i],'_')[[1]][1]
  curation <- str_replace(str_split(files_[i],'_')[[1]][5],"jaccard","")
  
  if(str_count(files_[i],"Rahman2015")>0){
    annotation_names_col <-TRUE
  } else {
    annotation_names_col <- FALSE
  }
  out2 <- pheatmap(jaccard_sim_, show_colames=FALSE,show_rownames=FALSE,labels_col =rep(' ',nrow(jaccard_sim_)),
                   #annotation_legend = T,annotation_names_col = T, 
                   main=paste0(model,"_",data,"_",curation), #fontsize_row = 3,
                   clustering_method = 'complete',cex=1,# color = magma(begin = 0, end = 1, n = 10),
                   #clustering_distance_rows = 'binary', clustering_distance_cols = 'binary',
                   cluster_cols=T, cluster_rows=T, scale="none",#cutree_rows = 2,
                   color= colorRampPalette(brewer.pal(n = 8, name = "Blues"))(length(mat_breaks)),#inferno(length(mat_breaks)),blues(length(mat_breaks)),
                   breaks = mat_breaks,border_color=F,annotation_colors = annotation_colors,#legend = F,
                   annotation_col=meta_samples,fontsize= 11, treeheight_row = 0, treeheight_col = 50
  )
  filename <- paste0("Figure/QC/" ,paste0(model,"_",data,"_",curation),"_Similarity_SamplLevel.png")
  #filename_2 <- paste0("Figure/QC/" ,paste0(model,"_",data,"_",curation),"_Similarity_SamplLevel.pdf")
  
  save_pheatmap_png(out2,filename)
  save_pheatmap_pdf(out2,filename_2)
  plot_list[[i]] <- out2[[4]]     ##to save each plot into a list. note the [[4]]
}
#wplot_save_this()
save_pheatmap_pdf <- function(x, filename, width=10, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

#save_pheatmap_pdf(plot_list[[12]], "test.pdf")
#save_pheatmap_png(plot_list[[12]],"filename.png",)

#plot(as.grob(out2))
library(sf)
library(cowplot)
#draw_plot(plot(as.grob(out2)), x = 2, y = 0, scale = 1)
# plot(out2$tree_col)
# cutree(as.hclust(out2$colDendrogram), 1)
# out2$tree_col %>%
#   as.dendrogram() %>%
#   plot(horiz = F)
# 
# # Create user-defined function, which extracts legends from ggplots
# extract_legend <- function(my_ggp) {
#   step1 <- ggplot_gtable(ggplot_build(my_ggp))
#   step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
#   step3 <- step1$grobs[[step2]]
#   return(step3)
# }
# 
# 
# shared_legend <- extract_legend(out2)
# z <- grid.arrange(arrangeGrob(plot_list[[1]] , plot_list[[3]],plot_list[[2]]  , ncol = 3),
#                    nrow = 2, heights = c(10, 10))
# z <- do.call(grid.arrange(shared_legend),plot_list)
# plot(z)
# ggsave("Figure/QC/g.png",z, units="in", width=25, height=15, dpi=150)


# 
# library(lemon)
# ## Heatmaps for all glioma classes excluding mixed-glioma AND GBM_IDH_wt
# dir.create(file.path('Figs/Revisiting_TCGA/', 'AST_vs_ODG'))
# for(i in 1:length(files)) {
#   # models_rxns <-read.csv(str_c('models/',files[i]),sep=',')
#   # rowname= models_rxns[,length(models_rxns)]
#   # rownames(models_rxns) = substr(rowname,1,nchar(rowname)-4)
#   # models_rxns = models_rxns[,1:length(models_rxns)-1]
#   # models_rxns = models_rxns[!rowSums(models_rxns)==0,]
#   # #models_rxns = models_rxns[!rowMeans(models_rxns)==1,]
#   # jaccard <- function(x, y) { return(length(intersect(x, y)) / length(union(x, y))) }
#   # jaccard_sim <- 1-as.matrix(dist(models_rxns, method = "binary") )
#   # jaccard_sim <- as.data.frame(jaccard_sim)
#   # #jaccard_sim$rownames <- rownames(models_rxns) 
#   # #colnames(jaccard_sim) <- substr(colnames(jaccard_sim) ,1,15)
#   # #rownames(jaccard_sim) <- substr(rownames(jaccard_sim) ,1,15)
#   jaccard_sim <-read_csv(str_c('models/',files[i]))
#   jaccard_sim <- as.data.frame(jaccard_sim)
#   rownames(jaccard_sim) <- jaccard_sim$Row 
#   rownames(jaccard_sim) <- str_replace(rownames(jaccard_sim),'.mat','')
#   colnames(jaccard_sim) <- str_replace(colnames(jaccard_sim),'.mat','')
#   jaccard_sim <- jaccard_sim[ , -which(names(jaccard_sim) %in% c("Row"))]
#   jaccard_sim <- 1-jaccard_sim
#   #
#   ## Selected intersected samples between the jaccard similarity and metadata
#   if (nchar(colnames(jaccard_sim)[1])==28) {
#     column <- "barcode"
#   } else {
#     column <- "sample"
#   }
#   meta_who2021 <- meta[meta$WHO_2021!="NA_",]
#   meta_who2021 <- meta_who2021[meta_who2021$WHO_2021!="ODG_AST",] # remove Mixed glioma
#   meta_who2021 <- meta_who2021[meta_who2021$WHO_2021!="GBM_IDH_wt",] # remove Mixed glioma
#   
#   intersected_samples <- intersect(rownames(jaccard_sim),meta_who2021[,column])
#   meta_samples <- meta_who2021[meta_who2021[,column] %in% intersected_samples,] 
#   
#   jaccard_sim_ <- jaccard_sim[rownames(jaccard_sim) %in% intersected_samples, colnames(jaccard_sim) %in% intersected_samples]
#   meta_samples<- meta_samples[order(rownames(jaccard_sim_)),]
#   rownames(meta_samples) <- rownames(jaccard_sim_)
#   meta_samples <- meta_samples[,c('paper_IDH.status','paper_X1p.19q.codeletion','WHO_2021')]#,'OS','Age',,'Study'
#   #mat_breaks <- quantile_breaks(jaccard_sim, n = 11)
#   #jaccard_sim$rownames <- NULL
#   data <- str_split(files[i],'_')[[1]][4]
#   model <- str_split(files[i],'_')[[1]][1]
#   curation <- str_replace(str_split(files[i],'_')[[1]][5],"jaccard","")
#   out2 <- pheatmap(jaccard_sim_, show_colames=FALSE,show_rownames=FALSE,labels_col =rep(' ',nrow(jaccard_sim_)),
#                    main=paste0(model,"_",data,"_",curation),
#                    clustering_method = 'complete',cex=1,# color = magma(begin = 0, end = 1, n = 10),
#                    #clustering_distance_rows = 'binary', clustering_distance_cols = 'binary',
#                    cluster_cols=T, cluster_rows=T, scale="none",#cutree_rows = 2,
#                    color= colorRampPalette(brewer.pal(n = 8, name = "Blues"))(length(mat_breaks)),#inferno(length(mat_breaks)),blues(length(mat_breaks)),
#                    breaks = mat_breaks,border_color=F,annotation_colors = annotation_colors,
#                    annotation_col=meta_samples,fontsize= 12, treeheight_row = 0, treeheight_col = 50
#   )
#   filename <- paste0("Figs/Revisiting_TCGA/AST_vs_ODG/" ,paste0(model,"_",data,"_",curation),"_Similarity_SamplLevel.png")
#   save_pheatmap_png(out2,filename)
#   #plot_list[[i]] <- out2[[4]]     ##to save each plot into a list. note the [[4]]
# }
# library(dendextend)
# 
# out2[[2]] %>%
#   as.dendrogram() %>%  
#   plot(horiz=F, lwd=2)
# # for (i in 1:length(files)) {
# #   data <- str_split(files[i],'_')[[1]][4]
# #   model <- str_split(files[i],'_')[[1]][1]
# #   curation <- str_replace(str_split(files[i],'_')[[1]][5],"jaccard","")
# #   png(filename=paste0("Figs/Revisiting_TCGA/" ,paste0(model,"_",data,"_",curation),"_Similarity_SamplLevel.png"), units="in", width=12, height=7, res=300)
# #   #print(plot_list[[i]])
# #   as.ggplot(plot_list[[i]] )
# #   dev.off()  
# # }
# # p1 <- as.ggplot(plot_list[[1]] )+ theme(legend.position="none")
# # 
# # data <- str_split(files[1],'_')[[1]][4]
# # model <- str_split(files[1],'_')[[1]][1]
# # curation <- str_replace(str_split(files[1],'_')[[1]][5],"jaccard","")
# # png(filename=paste0("Figs/Revisiting_TCGA/" ,paste0(model,"_",data,"_",curation),"_Similarity_SamplLevel.png"), units="in", width=12, height=7, res=300)
# # p1
# # dev.off() 
# 
# ##### Calculate the WHO_2021 matching
# df_keep <- df
# 
# df$rownames <- substr(df$rownames,1,15)
# df$colnames <- substr(df$colnames,1,15)
# df[df$rownames!=df$colnames,] -> df # remove identical samples
# # add the diagnosis classes
# meta$rownames <- meta$sample 
# meta$WHO_2021_x <- meta$WHO_2021
# df <- left_join(df,meta[,c('rownames','WHO_2021_x')])
# meta$colnames <- meta$sample
# meta$WHO_2021_y <- meta$WHO_2021
# df <- left_join(df,meta[,c('colnames','WHO_2021_y')])
# df <- distinct(df)
# df$WHO_2021_matching <- "Not_Matching"
# pool_ <- df$WHO_2021_x == df$WHO_2021_y
# df$WHO_2021_matching[pool_] <- "Matching"
# 
# df_nona <- df[ df$WHO_2021_x!="NA_",]
# df_nona <- df_nona[ df_nona$WHO_2021_y!="NA_",]
# 
# 
# #df_nona <- df_nona[! is.na(df_nona$WHO_2021_y),]
# 
# ##
# ## Boxplot of the relationship between matching and similarity
# 
# p <- ggplot(df_nona, aes(x=model, y=similarity)) + 
#   #geom_violin(aes(fill=as.character(WHO_2021_matching)),position = position_dodge(0.9) ) + 
#   geom_violin(aes(color = WHO_2021_matching), trim = FALSE,position = position_dodge(0.9) ) +
#   geom_boxplot(aes(color = WHO_2021_matching), width = 0.15,position = position_dodge(0.9)) +
#   scale_color_manual(values = c( "#00AFBB","#E7B800"))+
#   facet_grid(WHO_2021_x ~ WHO_2021_y,scales="free") + # # data ~ model
#   theme(axis.text.x = element_text(angle = 90,size = 13),
#         axis.text.y = element_text(angle = 0,size = 12)
#         #plot. = element_text(angle = 0,size = 14,face = "bold")
#         )
#   #stat_summary(fun = "mean", geom = "line")
# 
# png(filename="Figs/Revisiting_TCGA/Similarity_Matching_Boxplot.png", units="in", width=16, height=12, res=600)
# p
# dev.off() 
# 
# min(df$similarity)
# max(df$similarity)
# 
# ## Difference of the similarity median of matched to mismatch in s
# 
# df%>%group_by(data ,model,WHO_2021_matching) %>% 
#   summarise(Median=median(similarity),Mean=mean(similarity),  Std=sd(similarity)) -> df_median
#   
# df_median %>%group_by(data ,model) %>%
#   mutate(Median_Diff = Median[WHO_2021_matching=="Matching"] - Median[WHO_2021_matching=="Not_Matching" ])-> df_median
# 
# df_median[,4:7] <- signif(df_median[,4:7], digits = 3)
# 
# ##
# x <- df[df$data=='RSEM',]
# y <- df[df$data!='RSEM',]
# 
# # Adding the class names for y
# meta$rownames <- meta$barcode
# y <- left_join(y,meta[,c('rownames','TYPE')])
# meta$colnames <- meta$barcode
# meta$TYPE2 <- meta$TYPE
# y <- left_join(y,meta[,c('colnames','TYPE2')])
# 
# # Adding the class names for x
# meta = read.csv('data/TCGA_CGGA_metadata.csv')
# #meta$Cohort <- substr(meta$Study,1,4)
# meta$CGGA_ID = gsub('\\.','_',meta$CGGA_ID)
# meta$rownames <- meta$CGGA_ID
# x <- left_join(x,meta[,c('rownames','TYPE')])
# meta$colnames <- meta$CGGA_ID
# meta$TYPE2 <- meta$TYPE
# x <- left_join(x,meta[,c('colnames','TYPE2')])
# 
# 
# # append x and y after adding the class names
# df2 <- bind_rows(x,y)
# df2 <- df2[df2$TYPE!='UNK',] 
# df2 <- df2[df2$TYPE2!='UNK',]
# df2$TYPE_vs_TYPE <- paste(df2$TYPE,df2$TYPE2, sep='_vs_')
# unique(df2$TYPE_vs_TYPE )
# 
# df2 <- df2[df2$TYPE_vs_TYPE %in% c("AST_vs_GBM", "AST_vs_ODG" ,"ODG_vs_GBM"),]
# df2 <- df2[df2$rownames!=df2$colnames,]
# library(ggplot2)
# # Basic violin plot
# p <- ggplot(df2, aes(x=model, y=similarity,fill=data)) + 
#   geom_violin() + facet_grid(.~TYPE_vs_TYPE,scales="free") + 
#   stat_summary(fun = "mean", geom = "line")
# #+ geom_boxplot(width=0.1)
# p 
# p + stat_summary(fun.data="mean_sdl", mult=1, 
#                  geom="crossbar", width=0.2 )
# 
# 
# 
# df2 %>% group_by(data, model, TYPE_vs_TYPE)%>%
#   mutate(m = mean(similarity)) %>% distinct(data, model, TYPE_vs_TYPE,m)->df3
#   
# p <- ggplot(df3, aes(x=model, y=similarity,fill=data))+
#   geom_line(stat = "summary", fun = "mean")+
#   facet_grid(TYPE_vs_TYPE~.,scales="free")
#   
# p
# 
# # df2$rc <- paste(df2$rownames,df2$colnames, sep='_')
# # df2$cr <- paste(df2$colnames,df2$rownames, sep='_')
# # df2[df2$rc==df2$cr,] -> NULL
# # 
# # df2[paste0(df2[,c('rownames','colnames')]) == df2[,c('colnames','rownames')],]->n
# # n
# 
# 
# 
# ### Calculate accuracy
# df_nona$curation[df_nona$curation=='NoCuration'] <- 'NoMedium'
# df_nona%>% group_by(data,model,curation,rownames) %>%
#   arrange(desc(similarity),.by_group=TRUE) %>%
#   top_n(1,similarity)-> df_nona_top
# df_nona_top %>% group_by(data,model,curation) %>%
#   summarise(accuracy = n_distinct(rownames[WHO_2021_matching=='Matching'])/n_distinct(rownames)) ->df_nona_top_acc
# 
# 
# #files = list.files('models/',pattern = '_Model_Rxns')
# files = list.files('models/',pattern = 'statistics')
# files <- str_subset(files,'FPKM')
# files <- str_subset(files,'Subtypes',negate = TRUE)
# # Select 
# #patterns <- c('DrugMetabolismRemoved','NoMedium','CSF_Thiele2020'
# patterns <- c('NoCuration','CSF_Thiele2020')
# 
# Pattern = paste(patterns, collapse="|")
# files <- files[grepl(Pattern,files)]
# 
# df_rxns <- data.frame(matrix(ncol = 4, nrow = 0))
# colnames(df_rxns) <- c("data",'model',"curation","Median_reactions")
# df_rxns <- df_rxns %>%
#   mutate_if(is.logical, as.character)
# df_rxns$Median_reactions <- as.numeric(df_rxns$Median_reactions )
# 
# # Iterate for each model-data reactions and calculate the similarity
# for(i in 1:length(files)) {
#   stat_df <-read_csv(str_c('models/',files[i]))
#   stat_df$data <- str_split(files[i],'_')[[1]][4]
#   stat_df$model <- str_split(files[i],'_')[[1]][1]
#   stat_df$curation <- str_replace(str_split(files[i],'_')[[1]][5],"statistics.","")
#   stat_df$Median_reactions <- median(stat_df$Reactions,na.rm = T)
#   stat_df$Median_reactions <- as.numeric(stat_df$Median_reactions )
#   
#   stat_df <-stat_df[,colnames(df_rxns)]
#   df_rxns <- bind_rows(df_rxns,distinct(stat_df))
# }
# df_rxns$curation[df_rxns$curation=='NoCuration'] <- 'NoMedium'
# ## Merge model median reactions with accurary
# df_summary <- left_join(df_nona_top_acc,df_rxns)
# df_summary$curation[df_summary$curation=='NoMedium'] <- 'No Medium'
# 
# df_summary$data[df_summary$data=='Biolinks'] <- 'TCGABiolinks'
# df_summary$data[df_summary$data=='RSEM'] <- "Ceccarelli2016"
# df_summary$data[df_summary$data=='Rahman2015'] <- "Rahman2015"
# 
# df_summary$label <- ''
# df_summary$label[df_summary$data=='Rahman2015' & df_summary$curation=='CSF' & df_summary$model=='Recon3D'] <- "Recon3D + CSF\n+ Rahman et al 2015"
# 
# ggplot(df_summary,aes(x=accuracy,
#                          y=Median_reactions,color=model,alpha=curation,shape=data,label=label))+
#   #geom_tile(alpha=1)+
#   geom_point(size=4)+
#   ggtitle("Scatterplot between patients' model specificty and separability") +
#   geom_text_repel(box.padding = 0.7, max.overlaps = Inf)+
#   xlab("Model separability: the accuracy (%) of matching the glioma subtype\nclass to the nearest sample using euclidean similarity") +
#   ylab("Model specificity: the median number of reactions") +
#   theme_classic() +
#   scale_alpha_discrete(range = c( 1,0.5))+
#   scale_color_brewer(palette="Dark2",direction = -1)+
#   theme(axis.text.x = element_text(angle = 0,size = 12),
#         axis.text.y = element_text(angle = 0,size = 13))
