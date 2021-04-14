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
object <- readRDS("data/MCL_AIM_74_20210402_SCT_min_dist=0.5_spread=1.2.rds")
object@meta.data$SCT_snn_res.0.8 %<>% as.character %>% as.integer
idents.1 = sort(unique(object$SCT_snn_res.0.8))

Idents(object) = "SCT_snn_res.0.8"
cluster_markers = FindMarkers.UMI(object = object,ident.1 = idents.1[args],
                                  group.by = "SCT_snn_res.0.8",
                                  logfc.threshold = 0.1,
                                  only.pos = F,
                                  test.use = "MAST",
                                  latent.vars = "nFeature_SCT")
id = idents.1[args]
if(id < 10) id %<>% paste0("0",.)
write.csv(cluster_markers,file = paste0(path,id,"-markers_FC0.1_",idents.1[args],".csv"))