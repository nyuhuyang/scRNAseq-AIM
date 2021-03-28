# test
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","sctransform",
                   "harmony"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

test_df = data.frame(spread = c(rep(3:7/5,times = 10),rep(seq(1.5,16,by = 0.5), times= 2)),
                     min_dist = c(rep(1:10/10,each = 5),rep(3:4/10,each = 30)))
# Selina chose args = 8,9,12,18,48
print(spread <- test_df[args,"spread"])
print(min.dist <- test_df[args,"min_dist"])

load(file = "data/MCL_AIM_74_20210311_SCT.Rda")
DefaultAssay(object) = "SCT"
object@assays$SCT@scale.data = matrix(0,0,0)
cc_genes = unlist(cc.genes)
names(cc_genes) = NULL
more_VariableFeatures <- unique(c(VariableFeatures(object),cc_genes))
length(more_VariableFeatures)
VariableFeatures(object) <- more_VariableFeatures
length(VariableFeatures(object))

object %<>% ScaleData
npcs = 85
object %<>% RunPCA(verbose = T,npcs = npcs)

system.time(object %<>% RunHarmony.1(group.by = "orig.ident", dims.use = 1:npcs,
                                     theta = 2, plot_convergence = TRUE,
                                     nclust = 50, max.iter.cluster = 100))

object %<>% FindNeighbors(reduction = "harmony",dims = 1:npcs)
object %<>% FindClusters(resolution = 0.8)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))
object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs,min.dist = min.dist,spread = spread)

file.name = paste0("min_dist=",min.dist,"_spread=",spread)
g1 <- UMAPPlot.1(object,group.by = "cell.types", title = paste0("min.dist = ",min.dist,
                                                              ", spread = ",spread),
           do.print = F,do.return = T)

jpeg(paste0(save.path, "/", file.name,"_cell.types.jpeg"), units= "in",width=10, height=7,res=600)
print(g1)
dev.off()

Rshiny_path <- paste0("Rshiny/MCL_AIM_74_20210327_",file.name,"/")
samples <-  c("All_samples")
PrepareShiny(object, samples = samples, Rshiny_path = Rshiny_path,reduction = "umap",
             verbose = T)

object@assays$RNA = NULL
object@assays$SCT@scale.data = matrix(0,0,0)
saveRDS(object, file = paste0("data/MCL_AIM_74_20210327_SCT_",file.name,".rds"))