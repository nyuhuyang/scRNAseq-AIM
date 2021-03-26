# test
invisible(lapply(c("Seurat","dplyr","harmony"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

test_df = data.frame(spread = rep(3:7/5,times = 10),
                     min_dist = rep(1:10/10,each = 5))
print(spread <- test_df[args,"spread"])
print(min.dist <- test_df[args,"min_dist"])

load(file = "data/MCL_AIM_74_20210311_SCT.Rda")
object@assays$SCT@scale.data = matrix(0,0,0)
cc_genes = unlist(cc.genes)
names(cc_genes) = NULL
more_VariableFeatures <- unique(c(VariableFeatures(object),cc_genes))
length(more_VariableFeatures)
VariableFeatures(object) <- more_VariableFeatures
length(VariableFeatures(object))
set.seed(101)

npcs = 85
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs,min.dist = min.dist,spread = spread)

file.name = paste0("min_dist=",min.dist,"_spread=",spread)
g1 <- UMAPPlot.1(object,group.by = "cell.types", title = paste0("min.dist = ",min.dist,
                                                              ", spread = ",spread),
           do.print = F,do.return = T)
g2 <- UMAPPlot.1(object,group.by = "SCT_snn_res.0.8", title = paste("min.dist =",min.dist,
                                                                    ", spread = ",spread),
                do.print = F,do.return = T)

jpeg(paste0(save.path, "/", file.name,"_cell.types.jpeg"), units= "in",width=10, height=7,res=600)
print(g1)
dev.off()

jpeg(paste0(save.path, "/", file.name,"_res.0.8.jpeg"), units= "in",width=10, height=7,res=600)
print(g2)
dev.off()