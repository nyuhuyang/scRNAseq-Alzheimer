invisible(lapply(c("Seurat","dplyr","ggplot2","scater","magrittr","pbapply",
                   "cowplot","data.table"), function(x) {
                       suppressPackageStartupMessages(library(x,character.only = T))
                   }))
DEGs <- readxl::read_excel("output/20230311/2022-07-06- in integrated_snn_res.0.2_Microglia 0_1_2_3_4_5_6_7_8_9_10_11_13_14_15_16_filter.xlsx", sheet = "all")


# generate expression txt file for GSEA analysis
#' @param df FindAllMarker results. data.frame
#' @param k an integer for the number of folds. createFolds argment
#' @param do.return TRUE/FALSE
#' @param continuous.label NULL/continuous label #http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html?_Phenotype_Labels
#' @example PrepareGSEA(object, k = 50, continuous.label = major_cells)
PreparePrerankGSEA <- function(df, names_from = "group",values_from = "avg_log2FC",gene = "genes",
                             continuous.label = NULL,file.name = NULL,save.path = NULL,Mouse2Human = c("toupper","biomaRt")[1],
                             do.return = FALSE,...){
    
    if(!any(class(df) %in% "data.frame")) stop("Incorrect input format, must be data.frame")
    set.seed(201)

    if(is.null(save.path)) save.path <- file.path("output",gsub("-","",Sys.Date()))
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
    if(is.null(file.name)) file.name = "GSEA_prerank_"
    
    df %>% select(c(gene,names_from,values_from)) %>% split(f = .[,names_from]) %>% 
        lapply(function(x){
        file_name <- file.path(save.path, paste0(file.name,x[1,2],".rnk"))
        x %<>% select(c(gene,values_from)) %>% arrange(desc(!!rlang::sym(values_from)))
        data.table::fwrite(x, file = file_name,
                           sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)
        x
    }) -> DEGs_list
    
    if(do.return) return(DEGs_list)
    
}
PreparePrerankGSEA(DEGs,names_from = "cluster",values_from = "avg_log2FC",gene = "gene")