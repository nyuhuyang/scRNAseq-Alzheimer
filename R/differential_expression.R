invisible(lapply(c("dplyr","magrittr","tidyr","openxlsx",#"Seurat","MAST","future",
                   "gplots"), function(x) {
                           suppressPackageStartupMessages(library(x,character.only = T))
                   }))

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#===========================================================
csv_names = list.files("output/20220621",pattern = "integrated_snn_res.0.2")
deg_list <- pbapply::pblapply(csv_names, function(csv){
        tmp <- read.csv(paste0("output/20220621/",csv),row.names = 1)
        tmp %<>% arrange(desc(avg_log2FC))
        tmp$gene = rownames(tmp)
        tmp
})

deg = bind_rows(deg_list)
deg1 = filter(deg, avg_log2FC > 0) %>% group_by(cluster) %>% top_n(50, avg_log2FC)
deg2 = filter(deg) %>% group_by(cluster) %>% top_n(100, abs(avg_log2FC))
deg_list = list(deg1,deg2,deg)
names(deg_list) = c("postive","positve_negative","all")
write.xlsx(deg_list, file = "output/20220621/DEGs_integrated_snn_res.0.2.xlsx",
           colNames = TRUE, borders = "surrounding")

#===================
#1 st
deg1 <- readxl::read_excel("output/20220706/1st analysis/2022-07-06- in integrated_snn_res.0.2_Microglia 0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16 .xlsx")
deg1 %<>% filter(p_val_adj < 0.05)
deg.p = filter(deg1, avg_log2FC > 0) %>% group_by(cluster) %>% top_n(50, avg_log2FC)
deg.both = filter(deg1) %>% group_by(cluster) %>% top_n(100, abs(avg_log2FC))
deg_list = list(deg.p,deg.both,deg1)
names(deg_list) = c("postive","positve_negative","all")
write.xlsx(deg_list, file = "output/20220706/1st analysis/2022-07-06- in integrated_snn_res.0.2_Microglia 0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16_filter.xlsx",
           colNames = TRUE, borders = "surrounding")

#===================
#2 nd 
deg2 <- readxl::read_excel("output/20220706/2nd analysis/2022-07-07- in integrated_snn_res.0.2_Microglia 0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16_19 .xlsx")
deg2 %<>% filter(p_val_adj < 0.05)
deg.p = filter(deg2, avg_log2FC > 0) %>% group_by(cluster) %>% top_n(50, avg_log2FC)
deg.both = filter(deg2) %>% group_by(cluster) %>% top_n(100, abs(avg_log2FC))
deg_list = list(deg.p,deg.both,deg2)
names(deg_list) = c("postive","positve_negative","all")
write.xlsx(deg_list, file = "output/20220706/2nd analysis/2022-07-07- in integrated_snn_res.0.2_Microglia 0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16_19_filter.xlsx",
           colNames = TRUE, borders = "surrounding")
