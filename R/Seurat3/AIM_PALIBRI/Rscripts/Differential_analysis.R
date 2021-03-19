########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
invisible(lapply(c("Seurat","dplyr","cowplot","fgsea",
                   "magrittr","data.table","future","ggplot2","tidyr"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
# change the current plan to access parallelization
plan("multiprocess", workers = 4)
plan()

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# Need 64GB
# load files

object = readRDS(file = "data/B_AIM_74_20210311_SCT.rds")
# Need 64GB
DefaultAssay(object) = "SCT"
Idents(object) = "X6clusters"
X6clusters = paste0("C",1:6)

cluster_markers = FindMarkers.UMI(object = object,ident.1 = X6clusters[args],
                                  group.by = "X6clusters",logfc.threshold = 0,
                                  only.pos = T,
                                  return.thresh = 1,
                                  test.use = "MAST",
                                  latent.vars = "nFeature_SCT")
write.csv(cluster_markers,file = paste0(path,"markers_",X6clusters[args],"csv"))