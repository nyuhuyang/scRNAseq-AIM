invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save.path <- paste0("output/",gsub("-","",Sys.Date()))
if(!dir.exists(save.path)) dir.create(save.path, recursive = T)

#======================================
object = readRDS(file = "data/MCL_AIM_93_20220519.rds")
opts = c("test_npcs_dist_spread","densmap")[1]
if(opts == "test_npcs_dist_spread"){
    object@meta.data = readRDS("output/MCL_AIM_93_20220519_metadata_v2.rds")
    test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                         spread = rep(c(0.6,1.0,1.4),times = 3),
                         npcs = rep(c(50,60,70,80,90,100),each = 9))
    test_df = bind_rows(list(test_df,test_df))
    test_df$nfeatures = rep(c(2000,3000),each =54)
    
    for(i in 73:90){
        nfeatures <- test_df[i,"nfeatures"]
        spread <- test_df[i,"spread"]
        min.dist <- test_df[i,"min_dist"]
        npcs <- test_df[i,"npcs"]
        save.path <- paste0("output/20220526/",nfeatures)
        umap.name = paste0("nfeatures",nfeatures,"npcs",npcs,"dist.",min.dist,"spread.",spread) %>% gsub("\\.","",.)
        title = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread)
        file.name = paste0("output/20220526/",nfeatures,"/umap_",title,".rds")
        umap  = readRDS(file.name)[[1]]
        umap@key = "UMAP_"
        colnames(umap@cell.embeddings) = c("UMAP_1","UMAP_2")
        object[["umap"]] = umap

        UMAPPlot.1(object, group.by = "X6cluster_AIM74",do.print = T,
                   raster=FALSE,no.legend = T,label = T, label.repel = T,alpha = 0.85,
                   title = title,file.name = paste0(title, ".jpeg"),
                   save.path = save.path)
        Progress(i,90)
    }

}


#########################
if(opts == "densmap"){
    object@meta.data = readRDS("output/MCL_AIM_93_20220519_metadata_v2.rds")
    
    test_df = data.frame(min_dist = rep(c(0.2,0.4,0.6),each = 3),
                         spread = rep(c(0.6,1.0,1.4),times = 3),
                         dens_lambda = rep(c(1,2,5),each = 9),
                         npcs = rep(c(70,80),each = 27))
    

    
    
    for(args in 1:54){
        min.dist <- test_df[args,"min_dist"]
        spread <- test_df[args,"spread"]
        dens_lambda <- test_df[args,"dens_lambda"]
        npcs <- test_df[args,"npcs"]
        file.name = paste0("npcs",npcs,"_dist.",min.dist,"_spread.",spread,"_dens_lambda.",dens_lambda)
        
        read.path <- "output/20220720/"
        save.path <- "output/20220721/"
        if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
        
        if(!file.exists(paste0(read.path,"umap_",file.name,".rds"))){
            print(args); 
    }
        reductions  = readRDS(paste0(read.path,"umap_",file.name,".rds"))
        object@reductions = reductions
        UMAPPlot.1(object, group.by = "X6cluster_AIM74",do.print = T,
                   raster=FALSE,no.legend = T,label = T, label.repel = T,alpha = 0.85,
                   title = file.name,file.name = paste0(file.name, ".jpeg"),
                   save.path = save.path)
        Progress(args,54)
    }
}
