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

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
# Need 64GB
# load files

object = readRDS(file = "data/B_AIM_74_20210311_SCT.rds")
# Need 64GB
DefaultAssay(object) = "SCT"

select_ident  = c("X6clusters","SCT_snn_res.0.8")[2]
# Cluster C6 ==== 
object@meta.data$X6clusters %<>% factor(levels = paste0("C",1:6))
Idents(object) = "X6clusters"
X6clusters = paste0("C",1:6)

# Cluster SCT_snn_res.0.8 ==== 
object@meta.data$SCT_snn_res.0.8 %<>% as.character %>% as.integer
SCT_snn_res = sort(unique(object$SCT_snn_res.0.8))

Idents(object) = "SCT_snn_res.0.8"
Idents(object) %<>% factor(levels = SCT_snn_res)

#(csv_list <- list.files("output/20210319",pattern = "markers_C",full.names = T))
(csv_list <- list.files("output/20210402",pattern = "markers_FC0.1",full.names = T))
deg_list <- lapply(csv_list, function(x){
    tmp = read.csv(x,row.names = 1)
    tmp$cluster = sub(".*_","",x) %>% sub("\\.csv","",.)
    tmp$gene = rownames(tmp)
    return(tmp)
})
gde.markers <- bind_rows(deg_list)
#write.csv(gde.markers, file = paste0(path,"DEG_B_AIM_74_20210311.csv"))
write.csv(gde.markers, file = paste0(path,"DEG_B_AIM_74_FC0.1_20210402.csv"))

gde.markers = read.csv(file = "output/20210319/DEG_B_AIM_74_20210311.csv")
#DoHeatmap.1======
(mito.genes <- grep(pattern = "^MT-", x = gde.markers$gene))
if(length(mito.genes)>0) gde.markers = gde.markers[-mito.genes,]
Top_n = 10
top = gde.markers %>% group_by(cluster) %>% top_n(Top_n, avg_logFC)

markers <- FilterGenes(object,c("CCND1","CD19","CD5","CDK4","RB1","BTK","SOX11"))

features = c(as.character(top$gene),
             tail(VariableFeatures(object = object), 2),
             markers)
length(features)
object %<>% ScaleData(features=features)
featuresNum <- make.unique(features, sep = ".")
object %<>% MakeUniqueGenes(features = features)
DoHeatmap.1(object =object, features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 0, group.bar = F, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0, label=F, cex.row= 2, legend.size = 0,width=10, height=6.5,
            title = paste("Top",Top_n,"DE genes of SCT_snn_res.0.8 in B and MCL cells"),
            file.name = paste0("Heatmap_SCT_snn_res.0.8.jpeg"))

DoHeatmap.1(object =object, features = featuresNum, Top_n = Top_n,
            do.print=T, angle = 45, group.bar = T, title.size = 20, no.legend = F,size=5,hjust = 0.5,
            group.bar.height = 0.05, label=T, cex.row= 2, legend.size = 0,width=10, height=6.5,
            title = paste("Top",Top_n,"DE genes of 6 sublcusters in B and MCL cells"),
            file.name = paste0("Heatmap_B_X6clusters_Label.jpeg"))
# fgsea
res = gde.markers
res = res[order(res["p_val_adj"]),]
head(res, 20)
(clusters <- unique(res$cluster))
hallmark <- fgsea::gmtPathways("../seurat_resources/msigdb/h.all.v6.2.symbols.gmt")
names(hallmark) = gsub("HALLMARK_","",names(hallmark))
names(hallmark) = gsub("\\_"," ",names(hallmark))
hallmark$`NF-kB signaling` =  read.delim("../scRNAseq-MCL/data/200222 NFKB pathway gene list.txt") %>%
    pull %>% as.character()
hallmark$`MYC TARGETS` = c(hallmark$`MYC TARGETS V1`,hallmark$`MYC TARGETS V2`)

# Now, run the fgsea algorithm with 1000 permutations:
set.seed(100)
fgseaRes = FgseaDotPlot(stats=res, pathways=hallmark,
                        padj = 0.25,pval = 0.05,
                        #order.yaxis.by = c("CD16 Monocytes","NES"),
                        decreasing = F,
                        #order.xaxis = c("CD14 Monocytes","CD16 Monocytes"),
                        Rowv = F,
                        Colv = F,
                        cols = c("blue","green","yellow","orange","chocolate1","red"),
                        size = " -log10(pval)", fill = "NES",
                        pathway.name = "Hallmark",rotate.x.text = T,
                        #title = paste0("Hallmark_B_X6clusters_noLabel.jpeg"),
                        font.xtickslab=14, font.main=12, font.ytickslab = 10,
                        font.legend = list(size = 12),font.label = list(size = 12),
                        do.return = T,
                        save.path = path,
                        do.print = T,
                        width = 5,height = 7,hjust = 0.75)
write.csv(fgseaRes, file = paste0(save.path,"FgseaDotPlot_FDR1_pval1.csv"))
# VolcanoPlots
gde = gde.markers[gde.markers$cluster %in% "CD16 Monocytes",]
g <- VolcanoPlots(gde, cut_off_value = 0.05, cut_off = "p_val_adj",
                  cut_off_logFC = 1,top = 15,
                  cols = c("#ba2832","#d2dae2","#2a71b2"),
                  cols.order = c('Upregulated','Stable','Downregulated'),
                  alpha=1, size=3,
                  legend.size = 12,legend.position="bottom",
                  color = "black", pch=21)
g = g + ggtitle("CD16 Monocytes vs CD14 Monocytes")
g = g + TitleCenter()
print(g)

jpeg(paste0(save.path,"VolcanoPlots_",file.name,".jpeg"), units="in", width=10, height=10,res=600)
print(g)
dev.off()
