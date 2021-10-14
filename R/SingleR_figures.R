# conda activate r4.0.3
library(Seurat)
library(magrittr)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
object = readRDS("data/Macrophages_5_20211112.rds")
pred1 = readRDS(file = "output/Macrophages_5_20211112_singleR_pred.rds")
pred2 = readRDS(file = "output/Macrophages_5_20211112_azimuth_pred.rds")


singlerDF = data.frame("label.blue_encode" = pred1$pruned.labels,
                       "label.human_brain" = pred2$pruned.labels,
                       row.names = rownames(pred))
table(rownames(pred1) == rownames(object@meta.data))
table(is.na(singlerDF$label.blue_encode))
table(is.na(singlerDF$label.human_brain))

singlerDF$label.blue_encode[is.na(singlerDF$label.blue_encode)]= "unknown"
singlerDF$label.human_brain[is.na(singlerDF$label.human_brain)]= "unknown"

##############################
# adjust cell label
##############################

##############################
# process color scheme
##############################
table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF)
lapply(c("label.blue_encode","label.human_brain"), function(g){
    UMAPPlot.1(object = object, label = T, label.repel = T,group.by = g,
               no.legend = T,cols = Singler.colors,
               pt.size = 0.1,label.size = 5,alpha = 0.85,
               do.print = T,do.return = F,
               title =paste("labeling by",g))
})

saveRDS(object, file = "data/Macrophages_5_20211112.rds")


# by barplot
cell_Freq <- table(object$label.singler) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% round(digits = 2) %>% scales::percent()
cols = ExtractMetaColor(object)
cell_Freq$cols = cols[cell_Freq$Var1]
cell_Freq = cell_Freq[order(cell_Freq$Var1),]

cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          ylab = "Cell Number",
          label = cell_Freq$Percent,
          lab.size = 3,
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of cell types in total 6 samples")+NoLegend()+
    theme(plot.title = element_text(hjust = 0.5,size=15),
          axis.text.x = element_text(vjust = 0.5))
dev.off()
