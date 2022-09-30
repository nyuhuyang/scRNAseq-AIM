########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
####################################
invisible(lapply(c("Seurat","dplyr","magrittr","tidyr","data.table",
                   "future","gplots"), function(x) {
        suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

#==============
opts = matrix(c())
csv_names = sapply(opts,function(x) paste)
csv_index = list.files("output/20220420",pattern = ".csv") %>% gsub("-.*","",.) %>% as.integer()
table(1:13 %in% csv_index)
csv_names = list.files("output/20220630",pattern = ".csv")
deg_list <- pbapply::pblapply(csv_names, function(csv){
    tmp <- read.csv(paste0("output/20220630/",csv),row.names = 1)
    tmp = tmp[tmp$p_val_adj < 0.05,]
    tmp$gene = rownames(tmp)
    #tmp %<>% group_by(cluster) %>% arrange(desc(avg_log2FC), .by_group = TRUE)
    tmp
})

deg = bind_rows(deg_list)
deg %<>% filter(p_val_adj < 0.05)
deg = deg %>% 
    group_by(cluster) %>% 
    top_n(2, avg_log2FC)

openxlsx::write.xlsx(deg, file = paste0(path,"SCT_snn_res.0.8_B_DEG.xlsx"),
           colNames = TRUE, borders = "surrounding")

