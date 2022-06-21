library(Seurat)
library(magrittr)
library(dplyr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_differential_expression.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

# Need 64GB ?
set.seed(101)
# SLURM_ARRAY_TASK_ID
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
if (length(slurm_arrayid)!=1)  stop("Exact one argument must be supplied!")
# coerce the value to an integer
args <- as.integer(as.character(slurm_arrayid))
print(paste0("slurm_arrayid=",args))

object <- readRDS("data/Macrophages_5_20211112.rds")
object@meta.data = readRDS(file = "shinyApp/Human_brain/meta_data.rds")

(step = c("resolutions")[1])

if(step == "resolutions"){# 32GB
    opts = c("0", "1", "2", "3","4", "5", "6", "7", 
             "8", "9", "10", "11","14", "15", "16")
    object %<>% subset(subset = label.human_brain.v2 == "Microglia" &
                           integrated_snn_res.0.2 %in% opts
                       )
    object$integrated_snn_res.0.2 %<>% droplevels
    #==========================
    Idents(object) = "integrated_snn_res.0.2"
    opt = opts[args]
    markers = FindMarkers_UMI(object, ident.1 = opt,
                              group.by = "integrated_snn_res.0.2",
                              assay = "SCT",
                              #min.pct = 0.01,
                              logfc.threshold = 0.1,
                                 only.pos = F#,
                                 #test.use = "MAST",
                                 #latent.vars = "nFeature_SCT"
                              )
    markers$cluster = opt
    num = opt
    if(args < 10) num = paste0("0",num)
    if(args < 100) num = paste0("0",num)

    arg = args
    if(args < 10) arg = paste0("0",arg)
    if(args < 100) arg = paste0("0",arg)

    write.csv(markers,paste0(path,arg,"_integrated_snn_res.0.2_",num, ".csv"))
}