library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

# ======== 2.1 =========== test with known markers==================
object <- readRDS("data/Macrophages_5_20211112.rds")

object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
Microglia_1 = Microglia & object$UMAP_1 > 0
Microglia_2 = Microglia & object$UMAP_1 <= 0 & object$UMAP_2 > 1
Microglia_3 = Microglia & object$UMAP_1 <= 0 & object$UMAP_2 <= 1

object@meta.data[Microglia_1,"label.human_brain"] = "Microglia 1"
object@meta.data[Microglia_2,"label.human_brain"] = "Microglia 2"
object@meta.data[Microglia_3,"label.human_brain"] = "Microglia 3"

#========= pathway ======
object <- readRDS("data/Macrophages_5_20211112.rds")
object@meta.data = readRDS(file = "shinyApp/Human_brain/meta_data.rds")

GeneSets = read.csv("doc/Gene lists for gene set enrichments.csv")
GeneSets %<>% df2list

httr::set_config(httr::config(ssl_verifypeer = FALSE))
human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
                          values = GeneSets[["Mouse_DAM"]] , mart = mouse,
                          attributesL = c("hgnc_symbol"), 
                          martL = human, uniqueRows=T)
rm = duplicated(genesV2[,1])
genesV2 = genesV2[!rm,]
GeneSets[["Mouse_DAM"]] = genesV2$HGNC.symbol
table(GeneSets[["Mouse_DAM"]]  %in% rownames(object))
names(GeneSets)[!names(GeneSets) %in% colnames(object@meta.data)]

for(i in which(!names(GeneSets) %in% colnames(object@meta.data))){
    print(names(GeneSets)[i])
    object %<>% AddModuleScore(features = GeneSets[i],
                               name = names(GeneSets)[i])
    colnames(object@meta.data) %<>% sub(paste0(names(GeneSets)[i],"1"),
                                        names(GeneSets)[i],.)
}

for(i in 11:13){
    object@meta.data[,paste0(names(GeneSets)[i],".pos")] = object@meta.data[,names(GeneSets)[i]] -min(object@meta.data[,names(GeneSets)[i]])
}

#============ cnv =================
meta.data = read.csv("shinyApp/Human_brain/cnv_meta_data.csv",row.names = 1)
meta.data$cnv_leiden %<>% factor(levels = sort(unique(meta.data$cnv_leiden)))
meta.data$cnv = as.factor(as.character(round(meta.data$cnv_score,digits = 6)))
for(i in 1:length(df_samples$sample.id)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(df_samples$sample[i])
    print(table(cells))
    meta.data[cells,"patient"] = as.character(df_samples$patient[i])
}
meta.data$patient %<>% factor(levels = c("C11","AD34","AD52","AD53"))
object@meta.data %<>% cbind(meta.data[,c("cnv_leiden","cnv_score","cnv","patient")])
object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings)

DefaultAssay(object) = "integrated"
object %<>% FindNeighbors(reduction = "cca.umap",dims = 1:2)
for(res in c(0.01,0.1,seq(from = 0.2, to = 2, by =0.2))){
    print(res)
    object %<>% FindClusters(resolution = res, graph.name = "integrated_snn")
}

DefaultAssay(object) = "SCT"
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
for(res in c(0.01,0.1,seq(from = 0.2, to = 2, by =0.2))){
    print(res)
    object %<>% FindClusters(resolution = res, graph.name = "SCT_snn")
}

object$orig.ident %<>% factor(levels = c("NBB_C11_AG","NBB_AD34_HIP",
                                         "NBB_AD52_HIP","NBB_AD52_SPG","NBB_AD53_HIP"))
object$condition %<>% factor(levels = c("CTRL (no mutation)","AD (no mutation)","AD (MAPK mutation)"))
s.genes <- cc.genes$s.genes
s.genes %<>% gsub("MLF1IP","CENPU",.)
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(object) = "SCT"
object <- CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
colnames(object@meta.data)[grep("Phase",colnames(object@meta.data))]="cell.cycle.phase"

UMAP_cnv = read.csv("shinyApp/Human_brain/cnv_umap.csv",row.names = 1) %>% as.matrix()
colnames(UMAP_cnv) = c("UMAP_1","UMAP_2")
rownames(UMAP_cnv) = rownames(meta.data)
object[["cnv.umap"]] <- CreateDimReducObject(embeddings = UMAP_cnv,
                                             key = "cnvUMAP_", assay = DefaultAssay(object))
saveRDS(object@meta.data, file = "shinyApp/Human_brain/meta_data.rds")


#======== rename ident =================
meta.data =  readRDS(file = "shinyApp/Human_brain/meta_data.rds")

df_annotation <- readxl::read_excel("doc/20220425_Annotation adjustments_LW.xlsx",
                                    sheet = "20220425_LW")
resolutions = paste0("SCT_snn_res.",c(1.2,2))
meta.data$label.human_brain.v2 = meta.data$label.human_brain

for(i in 1:length(resolutions)){
    keep = which(!is.na(pull(df_annotation[,resolutions[i]])))
    for(m in keep){
        cl = pull(df_annotation[m,resolutions[i]])
        orig.ident = pull(df_annotation[m,"orig.ident"]) %>% str_trim
        change_from = meta.data[,resolutions[i]] == cl &
                      meta.data[,"orig.ident"] == orig.ident &
                        meta.data[,"label.human_brain"] == "Microglia"
        change_to = pull(df_annotation[m,"label.human_brain.v2"])

        meta.data[change_from,"label.human_brain.v2"] = change_to
        print(paste ("nCell = ",length(which(change_from)),",",resolutions[i],"at",cl,"------->",change_to))
    }
}

saveRDS(meta.data, file = "shinyApp/Human_brain/meta_data.rds")

