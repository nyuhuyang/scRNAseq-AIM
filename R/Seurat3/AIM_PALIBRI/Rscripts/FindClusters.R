invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

resolutions = seq(0.4,1.2, by = 0.1)

object <- readRDS("data/B_AIM_74_20210311_SCT.rds")

object %<>% FindClusters(resolution = resolutions[args],method = "igraph",
                         algorithm = "Leiden")
UMAPPlot.1(object, group.by=paste0("SCT_snn_res.",resolutions[args]),
           pt.size = 0.3,label = T,
           label.repel = T,alpha = 0.9,
           do.return = F,
           no.legend = T,label.size = 4, repel = T, 
           title = paste("res =",resolutions[args]),
           file.name = paste0("SCT_snn_res=",resolutions[args],".jpeg"))
           do.print = T, save.path = paste0(path,"test_res"))
