invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

#======================================
object = readRDS(file = "data/MCL_AIM_93_20220519.rds")
meta.data = object@meta.data
meta.data$barcode  = rownames(meta.data)
meta.data$X6cluster_MCL61 = meta.data$label1.blue_encode
meta.data$X6cluster_MCL61.colors = meta.data$label1.blue_encode.colors
meta.data$X6cluster_AIM74 = meta.data$label1.blue_encode
meta.data$X6cluster_AIM74.colors = meta.data$label1.blue_encode.colors
X6cluster.colors = c('#40A635','#FE8205','#8861AC','#E83C2D',"#FB9A99","#FFFF99")#brewer.pal(n = 6,"Paired")

#=======MCL_61_20220331_metadata.rds===============
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
resolutions = c( 0.01, 0.1, 0.2, 0.5,0.8)
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}

UMAPPlot.1(object, do.print = T,raster=FALSE,cols = Singler.colors,label = T,label.repel = T,group.by = "SCT_snn_res.0.1",no.legend=T)
UMAPPlot.1(object, do.print = T,raster=FALSE,cols = Singler.colors,label = T,label.repel = T,group.by = "SCT_snn_res.0.2",no.legend=T)
UMAPPlot.1(object, do.print = T,raster=FALSE,cols = Singler.colors,label = T,label.repel = T,group.by = "SCT_snn_res.0.5",no.legend=T)


meta.data = object@meta.data
meta.data$X6cluster = meta.data$cell.types
meta.data[meta.data$SCT_snn_res.0.1 %in% c("0","1","11","12","19"),"X6cluster"] = "1"
meta.data[meta.data$SCT_snn_res.0.1 %in% "10","X6cluster"] = "2"
meta.data[meta.data$SCT_snn_res.0.1 %in% "13","X6cluster"] = "3"
meta.data[meta.data$SCT_snn_res.0.1 %in% "14","X6cluster"] = "4"
meta.data[meta.data$SCT_snn_res.0.1 %in% "9","X6cluster"] = "5"
meta.data[meta.data$SCT_snn_res.0.1 %in% "18","X6cluster"] = "6"
meta.data[meta.data$SCT_snn_res.0.1 %in% "16","X6cluster"] = "Pt02"
meta.data$X6cluster %<>% gsub("Monocytes","Monocytes:CD14+",.)
meta.data[meta.data$SCT_snn_res.0.2 %in% "6" &
              meta.data$X6cluster %in% "Monocytes:CD14+","X6cluster"] = "Monocytes:CD16+"

saveRDS(meta.data, file = "data/MCL_61_20220331_metadata.rds")


#=======MCL_AIM_93_20220519_metadata_v2.rds===============
meta_data1 = readRDS("../scRNAseq-MCL/data/MCL_61_20220331_metadata.rds")
for( i in 1:6) {
    cells = rownames(meta_data1)[meta_data1$X6cluster %in% as.character(i)]
    print(table(meta.data[cells,"label1.blue_encode"]))
    meta.data[cells,"X6cluster_MCL61"] = paste0("C",as.character(i))
    meta.data[cells,"X6cluster_MCL61.colors"] = X6cluster.colors[i]
    
}
# copy X6cluster from MCL_61_20220331
MCL <- readRDS("data/B_AIM_74_20210311_SCT.rds")
meta_data2 = MCL@meta.data
for( i in 1:6) {
    cells = rownames(meta_data2)[meta_data2$X6cluster %in% paste0("C",as.character(i))]
    cells = cells[cells %in% rownames(meta.data)]
    print(table(meta.data[cells,"label1.blue_encode"]))
    meta.data[cells,"X6cluster_AIM74"] = paste0("C",as.character(i))
    meta.data[cells,"X6cluster_AIM74.colors"] = X6cluster.colors[i]
    
}

table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
object$response %<>% factor(levels = c("Normal","Untreated","CR","PR","PD"))
object$treatment %<>% factor(levels = c("Normal","Untreated","PALIBR+Ibrutinib"))

meta.data$X6cluster_MCL61 %<>% factor(levels = unique(c("MCL","unknown",sort(unique(meta.data$X6cluster_MCL61)))))
meta.data$X6cluster_MCL61 %<>% factor(levels = unique(c("MCL","unknown",sort(unique(meta.data$X6cluster_MCL61)))))

table(meta.data$X6cluster_MCL61,meta.data$X6cluster_MCL61.colors)
meta.data$X6cluster_MCL61.colors %<>% gsub("#ff0000","#CAB2D6",.)#Erythrocytes
meta.data$X6cluster_AIM74.colors %<>% gsub("#ff0000","#CAB2D6",.)#Erythrocytes
meta.data$label1.blue_encode.colors %<>% gsub("#ff0000","#CAB2D6",.)#Erythrocytes

saveRDS(meta.data, "output/MCL_AIM_93_20220519_metadata_v2.rds")

#=======MCL_AIM_93_20220519_metadata_v3.rds===============
object = readRDS("data/MCL_AIM_93_20220519.rds")
object@meta.data = readRDS("output/MCL_AIM_93_20220519_metadata_v2.rds")

#umap
object[["umap"]] = NULL
file.name = paste0("output/20220526/3000/umap_npcs70_dist.0.2_spread.1.4.rds")
umap  = readRDS(file.name)[[1]]
umap@key = "UMAP_"
colnames(umap@cell.embeddings) = c("UMAP_1","UMAP_2")
object[["umap"]] <- umap

plots <- UMAPPlot(object,group.by = "X6cluster_AIM74")

jpeg(paste0(save.path, "/UMAP.jpg"), units="in", width=10, height=7,res=600)
print(plots+ rectangle(x_left = -3, x_right =15, y_bottom = -13, y_top = 4))
dev.off()
meta.data = object@meta.data
umap_cord = object[["umap"]]@cell.embeddings
umap_rectangle = umap_cord[,"UMAP_1"] > -3 & 
                 umap_cord[,"UMAP_1"] < 15 &
                 umap_cord[,"UMAP_2"] > -13 &
                 umap_cord[,"UMAP_2"] < 4
table(umap_rectangle)

# B_MCL
B_MCL = meta.data$X6cluster_AIM74 %in% c("B cells","MCL")
# cluster
meta.data$SCT_snn_res.0.8 %<>% as.character()
df = table(meta.data[umap_rectangle,"SCT_snn_res.0.8"]) %>% as.data.frame
df = df[order(df$Freq,decreasing = T),]
df1 = df[df$Freq > 50,]
df1$Var1 %<>% as.character()
df1$Var1
meta.data$X6cluster_AIM74_res.0.8 = meta.data$X6cluster_AIM74

for(cl in df1$Var1){
    change_from = B_MCL & umap_rectangle & meta.data$SCT_snn_res.0.8 == cl
    print(table(change_from))
    print(cl)
    
    meta.data[change_from,"X6cluster_AIM74_res.0.8"] = cl
}
table(meta.data$X6cluster_AIM74_res.0.8)
meta.data$X9cluster = meta.data$X6cluster_AIM74_res.0.8
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("36","15"),"X9cluster"] = "C3"
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("8"),"X9cluster"] = "C6"
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("9","7","37","23","38","3"),"X9cluster"] = "C3.5"
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("42","26"),"X9cluster"] = "C7"
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("2","24","34","58","43"),"X9cluster"] = "C1"
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("48","35","32","49","44","39","22"),"X9cluster"] = "C1.5"
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("40"),"X9cluster"] = "C5"
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("16","10"),"X9cluster"] = "C4"
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("18","11","59","6","61","0"),"X9cluster"] = "C4.5"
meta.data[B_MCL & umap_rectangle & meta.data$X6cluster_AIM74_res.0.8 %in% 
              c("51"),"X9cluster"] = "C7"
table(meta.data$X9cluster)
meta.data1 = meta.data[!duplicated(meta.data$X6cluster_AIM74),c("X6cluster_AIM74","X6cluster_AIM74.colors")]
meta.data$X9cluster.colors = meta.data$X9cluster
meta.data$X9cluster.colors %<>% plyr::mapvalues(from = meta.data1$X6cluster_AIM74,
                                               to = meta.data1$X6cluster_AIM74.color)
meta.data$X9cluster.colors %<>% gsub("C1.5","#9F941D",.)
meta.data$X9cluster.colors %<>% gsub("C3.5","#B84F6D",.)
meta.data$X9cluster.colors %<>% gsub("C4.5","#F26B63",.)
meta.data$X9cluster.colors %<>% gsub("C6","#F49E63",.)
meta.data$X9cluster.colors %<>% gsub("C7","#76776D",.)
table(meta.data$X9cluster.colors)
saveRDS(meta.data, "output/MCL_AIM_93_20220519_metadata_v3.rds")
