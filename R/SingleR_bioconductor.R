#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)
library(data.table)
library(Matrix)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ======= load single cell dataset =================
object <- readRDS("data/Macrophages_5_20211112.rds")

sce <- SingleCellExperiment(list(logcounts=object[["RNA"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ====== load human-multiple-cortical-areas-smart-seq =============
counts <- fread("../seurat_resources/azimuth/human_primary_motorcortex/matrix.csv",
                header=TRUE,sep=",",quote="",
                stringsAsFactors=FALSE)
meta.data <- fread("../seurat_resources/azimuth/human_primary_motorcortex/metadata.csv",
                   header=TRUE,sep=",",quote="",
                   stringsAsFactors=FALSE)

barcode = pull(counts[,1])
counts = t(counts[,-1] )
colnames(counts) =  barcode
counts_sparse <- Matrix(counts, sparse=TRUE)
RowMeans = rowMeans2(counts_sparse)
table(RowSum>0)
counts_sparse = counts_sparse[RowSum>0,]
libsizes <- colSums(counts_sparse)
size.factors <- libsizes/mean(libsizes)
Logcounts <- log1p(t(t(counts_sparse)/size.factors))

meta.data %<>% tibble::column_to_rownames(var = "sample_name")

human_brain <- SingleCellExperiment(list(logcounts=Logcounts),
                                    colData=DataFrame(meta.data))
human_brain = subset(human_brain, subclass_label !="")
head(rownames(human_brain))
saveRDS(object = human_brain, file = "../seurat_resources/azimuth/human_primary_motorcortex/human_brain.rds")
#human_brain = readRDS(file = "../seurat_resources/azimuth/human_primary_motorcortex/human_brain.rds")
# ====== load blue_encode reference =============
blue_encode <- BlueprintEncodeData()
head(rownames(blue_encode))
common <- Reduce(intersect, list(rownames(sce),
                                 rownames(blue_encode)
))
length(common)
table(blue_encode$label.fine)
system.time(trained <- trainSingleR(ref = blue_encode[common,],
                                    labels=blue_encode$label.fine))
system.time(pred <- classifySingleR(sce[common,], trained))
saveRDS(object = pred, file = "output/Macrophages_5_20211112_singleR_pred.rds")

#=========================
common <- Reduce(intersect, list(rownames(sce),
                                 rownames(human_brain)
))
length(common)
table(human_brain$subclass_label)
system.time(trained <- trainSingleR(ref = human_brain[common,],
                                    labels=human_brain$subclass_label))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/Macrophages_5_20211112_azimuth_pred.rds")
