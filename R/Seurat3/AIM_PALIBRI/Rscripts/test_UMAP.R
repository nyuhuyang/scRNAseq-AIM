# test
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","sctransform",
                   "harmony"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

test_df = data.frame(min_dist = c(rep(1:10/10,each = 5),rep(3:4/10,each = 30)),
                     spread = c(rep(3:7/5,times = 10),rep(seq(1.5,16,by = 0.5), times= 2)))

# Selina chose args = 18,24,48
print(spread <- test_df[args,"spread"])
print(min.dist <- test_df[args,"min_dist"])

load(file = "data/MCL_AIM_74_20210311.Rda")
DefaultAssay(object) = "SCT"
add_cc_genes = FALSE
if(add_cc_genes){
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
}
length(VariableFeatures(object))


object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs,min.dist = min.dist,spread = spread)

file.name = paste0("min_dist=",min.dist,"_spread=",spread)
g1 <- UMAPPlot.1(object,group.by = "SCT_snn_res.0.8", title = paste0("min.dist = ",min.dist,
                                                              ", spread = ",spread),
           do.print = F,do.return = T)

jpeg(paste0(save.path, "/", file.name,"_res.0.8.jpeg"), units= "in",width=10, height=7,res=600)
print(g1)
dev.off()

Rshiny_path <- paste0("Rshiny/MCL_AIM_74_20210327_",file.name,"_by_samples/")
samples <- c("N01","N02","N03","N04","PtU01","PtU02","PtU03","PtU04",
             "Pt2_30Pd","Pt10_LN2Pd","Pt11_LN1","Pt11_1","Pt11_14","Pt11_28",
             "Pt11_31","Pt13_BM1","Pt13_1a","Pt13_1b","Pt16_3Pd",
             "Pt17_LN1","Pt17_2","Pt17_7","Pt17_12","Pt17_31",
             "Pt19_BM2Pd","Pt20_1","Pt25_SB1","Pt25_1","Pt25_1_8","Pt25_24",
             "Pt25_25Pd","Pt25_AMB25Pd","Pt27_1","Pt27_1_8","Pt27_12",
             "Pt28_LN1","Pt28_1","Pt28_4","Pt28_28","PtB13_Ibp","PtB13_Ib1",
             "PtB13_IbR","AIM2_1","AIM2_56","AIM3_1","AIM3_56","AIM6_1",
             "AIM6_56","AIM6_140","AIM6_190","AIM10_1","AIM10_85",
             "AIM10_173","AIM13_BM1","AIM13_1","AIM13_12","AIM13_68",
             "AIM13_84","AIM13_140","AIM17_BM1","AIM17_1","AIM17_4",
             "AIM17_16","AIM17_40","AIM17_EOT","AIM18_1","AIM18_85",
             "AIM18_148","AIM24_BM1","AIM24_1","AIM24_4","AIM24_16",
             "AIM24_28","AIM24_59")
PrepareShiny(object, samples = samples, Rshiny_path = Rshiny_path,reduction = "umap",
             verbose = T)

Rshiny_path <- paste0("Rshiny/MCL_AIM_74_20210327_",file.name,"_all/")
samples <- c("All_samples")
PrepareShiny(object, samples = samples, Rshiny_path = Rshiny_path,reduction = "umap",
             verbose = T)

object@assays$RNA = NULL
object@assays$SCT@scale.data = matrix(0,0,0)
saveRDS(object, file = paste0("data/MCL_AIM_74_20210402_SCT_",file.name,".rds"))

PrepareShiny(object, samples = samples, Rshiny_path = Rshiny_path,
             reduction = "umap",split.by = "SCT_snn_res.0.8",
             verbose = T)
