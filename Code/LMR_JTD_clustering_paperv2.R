rm(list=ls()) ### clear variables at the start, so that old variables aren't carried over

rm(list=ls()) ### clear variables at the start, so that old variables aren't carried over


#Updated - 4/11/19

#Below is a sample set of code for the single cell lung organoid data on my (Sophia) computer. Change for your needs.
#Important: this code cannot be run straight through - there will be times when you need to make user decisions before preceding to the next step
#Important 2: Genes are different between mice and humans. For human cells, all letters are capatilized. Mouse genes have a mix of capital and lower case. R is case sensitive. Make sure you choose correctly. For example, human mitochondrial genes begin with "MT-" while mouse mitochondrial genes begin with "mt-"
#This is not all inclusive - refer to the Seurat website (Satija lab) for additional analyses/visualizations

#Set working directory
#This is the folder where all of your files are stored. This allows you to save and call files/objects by a single name rather than the entire file each time. i.e. "sample.Robj" vs "~/Desktop/Sample Folder/sample.Robj"
#setwd('G:/My Drive/March29')

#These are the programs necessary to run this code that are not included in base R
library(Seurat)
library(dplyr)

#Input files
#Files can be inputted in two ways. You can either upload a folder directly from the 10x genomics read or you can input the cell-gene count table. I will go through the input for each method, but the lung data uses 10x

#Input files from 10x
#You need to open the folder that contains three files - "barcodes.tsv", "matrix.mtx", and "genes.tsv" (sometimes saved as features.tsv, if so, you need to change the name to genes.tsv)
#The files will come zipped (have a .gz format) - you need to open the the zipped files in order to call them here
#The data directory is the name of the FOLDER containing these three files

P1 <- Read10X(data.dir = "/nfs/turbo/mcomm-shealabngsdata/globus/AGC/4809/Sample1/filtered_feature_bc_matrix/") 
P2 <- Read10X(data.dir = "/nfs/turbo/mcomm-shealabngsdata/globus/AGC/4809/Sample2/filtered_feature_bc_matrix/")
P3 <- Read10X(data.dir = "/nfs/turbo/mcomm-shealabngsdata/globus/AGC/4831/Sample1/filtered_feature_bc_matrix/") 
P4 <- Read10X(data.dir = "/nfs/turbo/mcomm-shealabngsdata/globus/AGC/4831/Sample2/filtered_feature_bc_matrix/")










#Input files from count matrix
#In this case you are opening a single FILE. This code is commented out as the lung data is in a 10x folder
#control.data <- read.table("Control.txt", header = TRUE, row.names = 1, sep = "\t", as.is = TRUE)
#TGF.data <- read.table("Control.txt", header = TRUE, row.names = 1, sep = "\t", as.is = TRUE)






#Create Seurat objects
#This formats the data into an object that can be processed by the Seurat algorithm
#You are now able to assign labels to the data while vastly reducing the file size
#Create the object, only include genes that are present in at least 3 cells

P1 <- CreateSeuratObject(counts=P1, min.cells = 3,min.features=200)
P2 <- CreateSeuratObject(counts=P2, min.cells = 3,min.features=200)
P3 <- CreateSeuratObject(counts=P3, min.cells = 3,min.features=200)
P4 <- CreateSeuratObject(counts=P4, min.cells = 3,min.features=200)


#Assign an identity to the cells - this is only necessary when you are combining more than one file/Seurat object
P1@meta.data$sample <- "D7 OVA"
P2@meta.data$sample <- "D7 MOG"
P3@meta.data$sample <- "D9 OVA"
P4@meta.data$sample <- "D9 MOG"


P1@meta.data$Day <- "D7"
P2@meta.data$Day <- "D7"
P3@meta.data$Day <- "D9"
P4@meta.data$Day <- "D9"


P1@meta.data$Antigen <- "OVA"
P2@meta.data$Antigen <- "MOG"
P3@meta.data$Antigen <- "OVA"
P4@meta.data$Antigen <- "MOG"

#In some cases you will want to filter out "cells" that are actually debris. This can be done by removing those cells that have a low gene count (in this case below 200) or those with a high percentage of mitochondrial cells (> 0.05)
#The values and this process are up to the users discretion


P1[["percent.mt"]] <- PercentageFeatureSet(P1, pattern = "^mt-")
P2[["percent.mt"]] <- PercentageFeatureSet(P2, pattern = "^mt-")
P3[["percent.mt"]] <- PercentageFeatureSet(P3, pattern = "^mt-")
P4[["percent.mt"]] <- PercentageFeatureSet(P4, pattern = "^mt-")




VlnPlot(P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(P3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


P1 <- subset(P1, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)
P2 <- subset(P2, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)
P3 <- subset(P3, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)
P4 <- subset(P4, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)


TNBC <- merge(x = P1, y = c(P2,P3,P4))

saveRDS(TNBC,file="/home/jtdecker/Aaron_data_full.rds")

library(plyr)
library(ggplot2)

library(cowplot)
theme_set(theme_cowplot())



ifnb.list <- SplitObject(TNBC, split.by = "sample")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20,k.filter=25)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"

#immune.combined<-pbmc


s.genes <- str_to_sentence(cc.genes$s.genes)
g2m.genes <- str_to_sentence(cc.genes$g2m.genes)
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

immune.combined <- ScaleData(immune.combined, vars.to.regress=c("S.Score","G2M.Score","nUMI","percent.mt"), verbose = TRUE)
immune.combined<-FindVariableFeatures(immune.combined)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = TRUE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined<-FindClusters(immune.combined,resolution=.3)

TNBC<-immune.combined

immune.combined<-TNBC

saveRDS(TNBC,file="Aaron data with UMAP.rds")

TNBC<-readRDS("Aaron data with UMAP.rds")


DefaultAssay(TNBC)<-"RNA"
TNBC<-NormalizeData(TNBC)

TNBC<-FindClusters(TNBC,verbose=TRUE)



TNBC.markers<-FindAllMarkers(TNBC,verbose = TRUE,only.pos = TRUE)

# Cluster 0 is Monocytes. Cd14 high
# Cluster 1 is monocytes. Prdx1, plp, Spp1 high
# Cluster 2 is Cd8 T cells. Cd3, IL7r, Cd8.
# Cluster 3 is macrphages, high in complement
# Cluster 4 is neutrophils. S100a9, S100a8
# Cluster 5 is NK cells.  High Gzma, Klra proteins
# Cluster 6 are fibroblasts.  High Col3a1, Col1a1, Dcn
# Cluster 7 are DC (antigen presenting cells). High Cd209a, H2 genes. CD11b+ Itgam
# Cluster 8 are ActvatedvCd4 T cells, cd3, cd4, il2ra, ifng.  Probably TH1
# Cluster 9 are ActvatedvCd4 T cells, cd3, cd4, il2ra, Foxp3, IL17. Probably TH17 and Treg
# Cluster 10 are NK cells again. KLr genes, Eomes.
# Cluster 11 are B cells Cd79 and Ig s,
# Cluster 12 are also B cells? Maybe immature B cells.  Cd24, Cd74, Clec9a == cDCs CD103
# Cluster 13 is another monocyte subset. Makes HP, Ccr2, 
# Cluster 14 are neurons? Tmem150c Immature DC
# Cluster 15 are plasmablast. Makes Igkv "Tnfrsf13b", "Slamf7" CD20 CD19



new.cluster.ids<-c("Monocytes",
                   "Macrophage",
                   "CD8 T",
                   "Complement Macrophage",
                   "Neutrophil",
                   "NK Cell",
                   "Fibroblast",
                   "Dendritic Cell",
                   "CD4 T","TH17",
                   "Helper NK",
                   "Mature B Cell",
                   "Antigen-presenting B Cell",
                   "Chemoattractant Monocyte",
                   "Immature DC",
                   "Plasma Cells")
names(new.cluster.ids) <- levels(TNBC)
TNBC <- RenameIdents(TNBC, new.cluster.ids)


#Laila Code starts here ----
getwd()
setwd("/Users/lailarad/Documents/ProgEAE_scRNAseq/")
path = "/Users/lailarad/Documents/ProgEAE_scRNAseq/"
TNBC<-readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/Aaron data with UMAP.rds.crdownload")
library(ggprism)
library(Seurat)
library(Seurat,lib.loc = .libPaths()[2])
packageVersion("Seurat")
library(dplyr)
library(plyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(SeuratObject,lib.loc = .libPaths()[2])

Idents(TNBC)
DimPlot(TNBC, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 5, repel = T)
DimPlot(TNBC, reduction = "pca", label = TRUE, pt.size = 0.1, label.size = 5, repel = T)

DimPlot(TNBC, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T, group.by = "Day")
dev.off()
TNBC[[]]
TNBC@meta.data

new.cluster.ids<-c("Monocytes",
                   "Macrophages",
                   "CD8 T",
                   "Macrophages",
                   "Neutrophil",
                   "NK Cell",
                   "Fibroblast",
                   "Cd11b+ DC",
                   "CD4 T","CD4 T",
                   "Helper NK",
                   "B Cell",
                   "CD11b- CD103+ DCs",
                   "Monocytes",
                   "Il4i1+ DC",
                   "Plasmablast")


#names(new.cluster.ids) <- levels(TNBC)
Idents(TNBC) = "integrated_snn_res.0.3"
TNBC[[]]
TNBC <- RenameIdents(TNBC, '0' = 'Mono',
                     "1" = "Mac",
                     "2" = "CD8 T",
                     "3" = "Mac",
                     "4" = "Neut",
                     "5" = "NK",
                     "6" = "Fibro",
                     "7" = "CD11b+ DC",
                     "8" = "CD4 T",
                     "9" = "CD4 T",
                     "10" = "Helper NK",
                     "11" = "B Cell",
                     "12" = "CD103+ DC",
                     "13" = "Mono",
                     "14" = "Il4i1+ DC",
                     "15" = "Plasmablast")

TNBC <- RenameIdents(TNBC, '0' = 'Monocytes',
                     "1" = "Macrophage",
                     "2" = "CD8 T",
                     "3" = "Complement Macrophage",
                     "4" = "Neutrophil",
                     "5" = "NK Cell",
                     "6" = "Stromal Cell",
                     "7" = "CD11b+ DC",
                     "8" = "CD4 T",
                     "9" = "CD4 T",
                     "10" = "Helper NK",
                     "11" = "B Cell",
                     "12" = "CD103+ DC",
                     "13" = "Chemoattractant Monocyte",
                     "14" = "Il4i1+ DC",
                     "15" = "Plasmablast")


TNBC$new.cluster.ids = Idents(TNBC)
getwd()

# Set cell identity classes using SetIdent
DimPlot(TNBC, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T)
cells.use <- WhichCells(TNBC, expression=`Cd8a` > 0 |`Cd8b1` > 0 , idents = c("CD4 T"))
TNBC <- SetIdent(TNBC, cells = cells.use, value = 'CD8 T')
TNBC$new.cluster.ids = Idents(TNBC)
DimPlot(TNBC, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T)


saveRDS(TNBC,file="LMR_16renamedclusters.rds")


#TNBC<-readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/LMR_16renamedclusters.rds")

TNBC[[]]
library(ggpubr)
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

Idents(TNBC) = TNBC$sub_cluster
png(file = "UMAP_20clusters.png",
    units = "in",width = 11, height = 11, res = 400)
DimPlot(TNBC, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T, 
        group.by = "sub_cluster",
        cols = getPalette(20)) +
  theme_prism(base_family = "Arial", base_size = 18) + theme(legend.text = element_text(size = 28)) +
  labs( title = "UMAP of Cell Types") +
  NoLegend()
dev.off()


TNBC<-RunTSNE(TNBC, reduction = "pca", dims = 1:20)

png(file = "TSNE_clusters20.png",
    units = "in",width = 11, height = 11, res = 400)
DimPlot(TNBC, reduction = "tsne", label = T, pt.size = 0.5, label.size = 5, label.box = F,
        #label.color = "white",
        repel = T, 
        group.by = "sub_cluster",
        #split.by = "Antigen",
        cols = getPalette(20)) +
  theme_prism(base_family = "Arial", base_size = 18, base_fontface = "bold") + theme(legend.text = element_text(size = 28)) +
  labs( title = "TSNE of Cell Types") + NoLegend()
ggsave(file = "TSNE_clusters20.svg", width = 11, height = 11, pointsize = 12)
dev.off()


# Frequency plots
CellTypeTable = as.data.frame(proportions(table(TNBC$sub_cluster,TNBC@meta.data$sample), margin = 2))
CellTypeTable$Freq = CellTypeTable$Freq*100
levels(CellTypeTable$Var1)

#getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#display.brewer.all()

png(file = "CellProportions_20clusters.png",
    units = "in",width = 7, height = 6.5, res = 400)
ggplot(CellTypeTable,aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism(base_family = "Arial", base_size = 16)+
  labs(x = "Condition", y= "Percent of Cells") + 
  scale_fill_manual(values = getPalette(20), name = "Cell Type") + 
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  )
ggsave(file = "CellProportions_20clusters.svg", width = 7, height = 6.5, pointsize = 12)
dev.off()


png(file = "UMAP_BMES.png",
    units = "in",width = 11, height = 11, res = 400)
DimPlot(TNBC, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T, cols = getPalette(16)) +
theme_prism(base_family = "Arial", base_size = 24) + theme(legend.text = element_text(size = 28)) +
NoLegend()
dev.off()
#ggsave(file = "DimPlotUMAP.svg", width = 14, height = 10, pointsize = 12)

png(file = "UMAP_BMES_day.png",
    units = "in",width = 11, height = 11, res = 400)
DimPlot(TNBC, reduction = "umap", label = F, pt.size = 0.5, group.by = "Day") +
  theme_prism(base_family = "Arial", base_size = 24) + theme(legend.text = element_text(size = 28)) 
dev.off()

png(file = "UMAP_BMES_Antigen.png",
    units = "in",width = 11, height = 11, res = 400)
DimPlot(TNBC, reduction = "umap", label = F, pt.size = 0.5, group.by = "Antigen") +
  theme_prism(base_family = "Arial", base_size = 24) + theme(legend.text = element_text(size = 28)) 
dev.off()

# Frequency plots
CellTypeTable = as.data.frame(proportions(table(Idents(TNBC),TNBC@meta.data$sample), margin = 2))
CellTypeTable$Freq = CellTypeTable$Freq*100
levels(CellTypeTable$Var1)

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
#display.brewer.all()

png(file = "CellProportionsBMES.png",
    units = "in",width = 7, height = 5, res = 400)
ggplot(CellTypeTable,aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism(base_family = "Arial", base_size = 16)+
  labs(x = "Condition", y= "Percent of Cells") + 
  scale_fill_manual(values = getPalette(16), name = "Cell Type") + 
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
        )
dev.off()

CellTypeData = CellTypeTable
# load stringr library
library(stringr)

# Split name column into firstname and last name
CellTypeData[c('Day', 'Antigen')] <- str_split_fixed(CellTypeData$Var2, ' ', 2)

plotCellType = function(data, celltype){
ggplot(data[which(data$Var1 == celltype),], aes(x = Antigen, y= Freq, fill = Day)) + 
  geom_bar(position = position_dodge(),stat = 'identity') +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism(base_family = "Arial", base_size = 5)+
  labs(x = "Condition", y= "Percent of Cells", title = paste(celltype)) + 
  theme(legend.text = element_text(size = 5)
  )
}

myplots = lapply(unique(CellTypeData$Var1), FUN = plotCellType, data = CellTypeData)

ggarrange(plotlist = myplots, ncol=4, nrow=4, 
          common.legend = T, legend = "right")
dev.off()

plotCellTypePresentation = function(data,celltype){
  ggplot(data[which(data$Var1 == celltype),], aes(x = Antigen, y= Freq, fill = Day)) + 
    geom_bar(position = position_dodge(),stat = 'identity') +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    theme_prism(base_family = "Arial", base_size = 16)+
    labs(x = "Condition", y= "Percent of Cells", title = paste(celltype)) + 
    theme(legend.text = element_text(size = 16)
    )
}

png(file = "Bcell_prop_BMES.png",
    units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentation(CellTypeData, "B Cell")
dev.off()
plotCellTypePresentation(CellTypeData, "CD8 T")
plotCellTypePresentation(CellTypeData, "Plasmablast")

myplots = lapply(c("B Cell", "CD8 T","CD4 T", "Plasmablast"), FUN = plotCellTypePresentation, data = CellTypeData)

png(file = "ImmuneCell_prop_BMES.png",
    units = "in",width = 21, height = 5, res = 400)
ggarrange(plotlist = myplots, ncol=4, nrow=1, 
          common.legend = T, legend = "right")
dev.off()


CellTypeCounts =as.data.frame(table(Idents(TNBC),TNBC@meta.data$sample))

# Split name column into firstname and last name
CellTypeCounts[c('Day', 'Antigen')] <- str_split_fixed(CellTypeCounts$Var2, ' ', 2)

myplots = lapply(unique(CellTypeData$Var1), FUN = plotCellType, data = CellTypeCounts)

ggarrange(plotlist = myplots, ncol=4, nrow=4, 
          common.legend = T, legend = "right")
dev.off()


#Determining cluster IDs ----

get_conserved <- function(cluster){
  FindConservedMarkers(TNBC,
                       ident.1 = cluster,
                       grouping.var = "Antigen",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    #left_join(y = unique(annotations[, c("gene_name", "description")]),
    #          by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}


conserved_markers <- map_dfr(unique(new.cluster.ids), get_conserved)
conserved_markers <- map_dfr(0:15, get_conserved)

DefaultAssay(immune.combined.sct)

# Extract top 10 markers per cluster
top20 <- conserved_markers %>% 
  mutate(avg_fc = (OVA_avg_log2FC + MOG_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)

TNBC[[]]
plots <- VlnPlot(TNBC, features = c("Cd3e", "Cd4", "Cd8a","Cd8b1", "Cd19", "Ms4a1", "Sdc1", "Cd27", "Slamf7"), 
                 split.by = "sample",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

plots <- VlnPlot(TNBC, features = c("Cd19", "Ms4a1", "Sdc1", "Tnfrsf13b", "Slamf7", "Ifngr1", "Cd74"), 
                 split.by = "sample",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)



FeaturePlot(TNBC, 
            reduction = "umap", 
            features = c("Cd4", "Cd8a", "Cd8b1"), 
            #sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

plots <- VlnPlot(TNBC, features = c(
  "Mcpt1", "Fcer1a"
), 
group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

plots <- VlnPlot(TNBC, features = c("Ptprc", "Flt1", "Ly6a",
                                                   "Id3", "Fhl2", "Gng11", "Cd34", "Acvrl1",
                                                   "Itgb3", "Vcam1", "Thbd", "Cd55", "Vegfa", "Vegfd",
                                                   "Ece1", "Plec", "Ecscr","Tgfbr2", "Icam1"), split.by = "Antigen",
                 group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

plots <- VlnPlot(TNBC, features = c("Cd14", "Lyz1"), 
                 split.by = "sample",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

#DCs
plots <- VlnPlot(TNBC, features = c("Cd209a", "Btla", "Flt3", "Itgax", #cd11c
                                    "Clec9a", "Ly75", "Sirpa", "Bst2", "Itgae", #CD103
                                    "Itgam", #Cd11b
                                    "Cd80", "Cd86", "Il4i1"
                                    ), 
                 split.by = "sample",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

plots <- VlnPlot(TNBC, features = c("Cd209a", "Itgax", #cd11c
                                    
                                    "Itgam", #Cd11b
                                    "Il4i1"
), 
split.by = "sample",
group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

plots <- VlnPlot(TNBC, features = c( "Btla", "Flt3", 
                                   "Ly75", "Sirpa", "Bst2", 
                                    
                                    "Il4i1"
), 
split.by = "sample",
group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

plots <- VlnPlot(TNBC, features = c("Flt3","H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa",
                                    "Cd80", "Cd86", "Cd83"
                                    ), 
                 split.by = "sample",
                 group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)
Idents(TNBC) = "seurat_clusters"
DCs_types_markers <- FindConservedMarkers(TNBC, ident.1 = "14", ident.2 =  c("7", "12"),grouping.var = "Antigen",
                                 verbose = FALSE)

plots <- VlnPlot(TNBC, features = c(
                                    "Cd80", "Cd86", "Ccl2", "Ccr2"
), 
split.by = "sample",
group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

TNBC[[]]
TNBC$celltype.sample <- paste(TNBC$new.cluster.ids, TNBC$sample, sep = "_")
Idents(TNBC) <- "celltype.sample"
unique(Idents(TNBC))
library(MAST)
TNBC
Infmono.de.D7 <- FindMarkers(TNBC, ident.1 = "Chemoattractant Monocyte_D7 MOG", ident.2 = "Chemoattractant Monocyte_D7 OVA", 
                             verbose = FALSE, test.use = "wilcox")
head(Infmono.de.D7, n = 25)
Infmono.de.D7["Ccl2",]
Infmono.de.D7["Ccr2",]
Infmono.de.D7["Cd80",]
Infmono.de.D7["Cd86",]
write.csv(Infmono.de.D7,paste(path,"/EAE_sc_Infmono.de.D7_MogvOVA.csv", sep = ""), row.names = T)

Infmono.de.D9 <- FindMarkers(TNBC, ident.1 = "Chemoattractant Monocyte_D9 MOG", ident.2 = "Chemoattractant Monocyte_D9 OVA", 
                             verbose = FALSE, test.use = "wilcox")
Infmono.de.D9["Ccl2",]
Infmono.de.D9["Ccr2",]
Infmono.de.D9["Cd80",]
Infmono.de.D9["Cd86",]
write.csv(Infmono.de.D9,paste(path,"/EAE_sc_Infmono.de.D9_MogvOVA.csv", sep = ""), row.names = T)


mono.de.D7 <- FindMarkers(TNBC, ident.1 = "Monocytes_D7 MOG", ident.2 = "Monocytes_D7 OVA", 
                          verbose = FALSE, test.use = "wilcox")
head(mono.de.D7, n = 25)
mono.de.D7["Ccl2",]
mono.de.D7["Ccr2",]
mono.de.D7["Cd80",]
mono.de.D7["Cd86",]
write.csv(mono.de.D7,paste(path,"/EAE_sc_mono.de.D7_MogvOVA.csv", sep = ""), row.names = T)


mono.de.D9 <- FindMarkers(TNBC, ident.1 = "Monocytes_D9 MOG", ident.2 = "Monocytes_D9 OVA", 
                          verbose = FALSE, test.use = "wilcox")
mono.de.D9["Ccl2",]
mono.de.D9["Ccr2",]
mono.de.D9["Cd80",]
mono.de.D9["Cd86",]
write.csv(mono.de.D9,paste(path,"/EAE_sc_mono.de.D9_MogvOVA.csv", sep = ""), row.names = T)


mac.de.D7 <- FindMarkers(TNBC, ident.1 = "Macrophage_D7 MOG", ident.2 = "Macrophage_D7 OVA", 
                         verbose = FALSE, test.use = "wilcox")
mac.de.D7["Ccl2",]
mac.de.D7["Ccr2",]
mac.de.D7["Cd80",]
mac.de.D7["Cd86",]

write.csv(mac.de.D7,paste(path,"/EAE_sc_mac.de.D7_MogvOVA.csv", sep = ""), row.names = T)


mac.de.D9 <- FindMarkers(TNBC, ident.1 = "Macrophage_D9 MOG", ident.2 = "Macrophage_D9 OVA", 
                         verbose = FALSE, test.use = "wilcox")
mac.de.D9["Ccl2",]
mac.de.D9["Ccr2",]
mac.de.D9["Cd80",]
mac.de.D9["Cd86",]
write.csv(mac.de.D9,paste(path,"/EAE_sc_mac.de.D9_MogvOVA.csv", sep = ""), row.names = T)


compmac.de.D7 <- FindMarkers(TNBC, ident.1 = "Complement Macrophage_D7 MOG", ident.2 = "Complement Macrophage_D7 OVA", verbose = FALSE, test.use = "wilcox")
compmac.de.D7["Ccl2",]
compmac.de.D7["Ccr2",]
compmac.de.D7["Cd80",]
compmac.de.D7["Cd86",]
write.csv(compmac.de.D7,paste(path,"/EAE_sc_compmac.de.D7_MogvOVA.csv", sep = ""), row.names = T)


compmac.de.D9 <- FindMarkers(TNBC, ident.1 = "Complement Macrophage_D9 MOG", ident.2 = "Complement Macrophage_D9 OVA", 
                             verbose = FALSE, test.use = "wilcox")
compmac.de.D9["Ccl2",]
compmac.de.D9["Ccr2",]
compmac.de.D9["Cd80",]
compmac.de.D9["Cd86",]
write.csv(compmac.de.D9,paste(path,"/EAE_sc_compmac.de.D9_MogvOVA.csv", sep = ""), row.names = T)


CD11bDC.de.D7 <- FindMarkers(TNBC, ident.1 = "CD11b+ DC_D7 MOG", ident.2 = "CD11b+ DC_D7 OVA", 
                             verbose = FALSE, test.use = "wilcox")
CD11bDC.de.D7["Ccl2",]
CD11bDC.de.D7["Ccr2",]
CD11bDC.de.D7["Cd80",]
CD11bDC.de.D7["Cd86",]
write.csv(CD11bDC.de.D7,paste(path,"/EAE_sc_CD11bDC.de.D7_MogvOVA.csv", sep = ""), row.names = T)


CD11bDC.de.D9 <- FindMarkers(TNBC, ident.1 = "CD11b+ DC_D9 MOG", ident.2 = "CD11b+ DC_D9 OVA", 
                             verbose = FALSE, test.use = "wilcox")
CD11bDC.de.D9["Ccl2",]
CD11bDC.de.D9["Ccr2",]
CD11bDC.de.D9["Cd80",]
CD11bDC.de.D9["Cd86",]
write.csv(CD11bDC.de.D9,paste(path,"/EAE_sc_CD11bDC.de.D9_MogvOVA.csv", sep = ""), row.names = T)


CD103DC.de.D7 <- FindMarkers(TNBC, ident.1 = "CD103+ DC_D7 MOG", ident.2 = "CD103+ DC_D7 OVA", 
                             verbose = FALSE, test.use = "wilcox")
CD103DC.de.D7["Ccl2",]
CD103DC.de.D7["Ccr2",]
CD103DC.de.D7["Cd80",]
CD103DC.de.D7["Cd86",]
write.csv(CD103DC.de.D7,paste(path,"/EAE_sc_CD103DC.de.D7_MogvOVA.csv", sep = ""), row.names = T)


CD103DC.de.D9 <- FindMarkers(TNBC, ident.1 = "CD103+ DC_D9 MOG", ident.2 = "CD103+ DC_D9 OVA", 
                             verbose = FALSE, test.use = "wilcox")
CD103DC.de.D9["Ccl2",]
CD103DC.de.D9["Ccr2",]
CD103DC.de.D9["Cd80",]
CD103DC.de.D9["Cd86",]
write.csv(CD103DC.de.D9,paste(path,"/EAE_sc_CD103DC.de.D9_MogvOVA.csv", sep = ""), row.names = T)



Il4i1DC.de.D7 <- FindMarkers(TNBC, ident.1 = "Il4i1+ DC_D7 MOG", ident.2 = "Il4i1+ DC_D7 OVA",
                             verbose = FALSE, test.use = "wilcox")
Il4i1DC.de.D7["Ccl2",]
Il4i1DC.de.D7["Ccr2",]
Il4i1DC.de.D7["Cd80",]
Il4i1DC.de.D7["Cd86",]
write.csv(Il4i1DC.de.D7,paste(path,"/EAE_sc_Il4i1DC.de.D7_MogvOVA.csv", sep = ""), row.names = T)


Il4i1DC.de.D9 <- FindMarkers(TNBC, ident.1 = "Il4i1+ DC_D9 MOG", ident.2 = "Il4i1+ DC_D9 OVA",
                             verbose = FALSE, test.use = "wilcox")
Il4i1DC.de.D9["Ccl2",]
Il4i1DC.de.D9["Ccr2",]
Il4i1DC.de.D9["Cd80",]
Il4i1DC.de.D9["Cd86",]
Il4i1DC.de.D9["Chi3l1",]
write.csv(Il4i1DC.de.D9,paste(path,"/EAE_sc_Il4i1DC.de.D9_MogvOVA.csv", sep = ""), row.names = T)

Il4i1DC.de.D9 <- FindMarkers(TNBC, ident.1 = "Il4i1+ DC_D9 MOG", ident.2 = "Il4i1+ DC_D9 OVA",
                             verbose = FALSE, test.use = "wilcox")
Il4i1DC.de.D9["Ccl2",]
Il4i1DC.de.D9["Ccr2",]
Il4i1DC.de.D9["Cd80",]
Il4i1DC.de.D9["Cd86",]
Il4i1DC.de.D9["Chi3l1",]
write.csv(Il4i1DC.de.D9,paste(path,"/EAE_sc_Il4i1DC.de.D9_MogvOVA.csv", sep = ""), row.names = T)

unique(TNBC$new.cluster.ids)

CD4T.de.D7 <- FindMarkers(TNBC, ident.1 = "CD4 T_D7 MOG", ident.2 = "CD4 T_D7 OVA",
                             verbose = FALSE, test.use = "wilcox")

CD4T.de.D9 <- FindMarkers(TNBC, ident.1 = "CD4 T_D9 MOG", ident.2 = "CD4 T_D9 OVA",
                          verbose = FALSE, test.use = "wilcox")

CD8T.de.D7 <- FindMarkers(TNBC, ident.1 = "CD8 T_D7 MOG", ident.2 = "CD8 T_D7 OVA",
                          verbose = FALSE, test.use = "wilcox")

CD8T.de.D9 <- FindMarkers(TNBC, ident.1 = "CD8 T_D9 MOG", ident.2 = "CD8 T_D9 OVA",
                          verbose = FALSE, test.use = "wilcox")

Stromal.de.D7 <- FindMarkers(TNBC, ident.1 = "Stromal Cell_D7 MOG", ident.2 = "Stromal Cell_D7 OVA",
                          verbose = FALSE, test.use = "wilcox")

Stromal.de.D9 <- FindMarkers(TNBC, ident.1 = "Stromal Cell_D9 MOG", ident.2 = "Stromal Cell_D9 OVA",
                             verbose = FALSE, test.use = "wilcox")

Neut.de.D7 <- FindMarkers(TNBC, ident.1 = "Neutrophil_D7 MOG", ident.2 = "Neutrophil_D7 OVA",
                             verbose = FALSE, test.use = "wilcox")

Neut.de.D9 <- FindMarkers(TNBC, ident.1 = "Neutrophil_D9 MOG", ident.2 = "Neutrophil_D9 OVA",
                             verbose = FALSE, test.use = "wilcox")

NK.de.D7 <- FindMarkers(TNBC, ident.1 = "NK Cell_D7 MOG", ident.2 = "NK Cell_D7 OVA",
                          verbose = FALSE, test.use = "wilcox")

NK.de.D9 <- FindMarkers(TNBC, ident.1 = "NK Cell_D9 MOG", ident.2 = "NK Cell_D9 OVA",
                          verbose = FALSE, test.use = "wilcox")

HNK.de.D7 <- FindMarkers(TNBC, ident.1 = "Helper NK_D7 MOG", ident.2 = "Helper NK_D7 OVA",
                        verbose = FALSE, test.use = "wilcox")

HNK.de.D9 <- FindMarkers(TNBC, ident.1 = "Helper NK_D9 MOG", ident.2 = "Helper NK_D9 OVA",
                        verbose = FALSE, test.use = "wilcox")

B.de.D7 <- FindMarkers(TNBC, ident.1 = "B Cell_D7 MOG", ident.2 = "B Cell_D7 OVA",
                         verbose = FALSE, test.use = "wilcox")

B.de.D9 <- FindMarkers(TNBC, ident.1 = "B Cell_D9 MOG", ident.2 = "B Cell_D9 OVA",
                         verbose = FALSE, test.use = "wilcox")

Plasma.de.D7 <- FindMarkers(TNBC, ident.1 = "Plasmablast_D7 MOG", ident.2 = "Plasmablast_D7 OVA",
                       verbose = FALSE, test.use = "wilcox")

Plasma.de.D9 <- FindMarkers(TNBC, ident.1 = "Plasmablast_D9 MOG", ident.2 = "Plasmablast_D9 OVA",
                       verbose = FALSE, test.use = "wilcox")

Top.Plasma.de.D7 = Plasma.de.D7[which(Plasma.de.D7$p_val_adj < 0.05 | abs(Plasma.de.D7$avg_log2FC)> 0.25),]
Top.Plasma.de.D7$avg_log2FC = abs(Top.Plasma.de.D7$avg_log2FC)
Top.Plasma.de.D7 = Top.Plasma.de.D7[ order(Top.Plasma.de.D7$p_val_adj,-Top.Plasma.de.D7$avg_log2FC),]
Top5genes_cluster_D7 = data.frame(Plasmablast = rownames(Top.Plasma.de.D7[1:5,]))

D7_markers_list = list(Plasma.de.D7, B.de.D7, HNK.de.D7, NK.de.D7, Neut.de.D7, Stromal.de.D7,
                       CD8T.de.D7, CD4T.de.D7, Il4i1DC.de.D7, mac.de.D7, Infmono.de.D7, mono.de.D7,
                       CD103DC.de.D7, CD11bDC.de.D7, compmac.de.D7)

Top5genes_cluster_D7 <- data.frame(matrix(ncol = 15, nrow = 5))
for (i in 1:length(D7_markers_list)){
  Top = D7_markers_list[[i]]
  Top$avg_log2FC = abs(Top$avg_log2FC)
  Top = Top[which(Top$p_val_adj < 0.05 | Top$avg_log2FC> 0.25),]
  Top = Top[ order(Top$p_val_adj,-Top$avg_log2FC),]
  Top = Top[!grepl("^Gm|Rik$|^Rps|^Rpl|^mt", rownames(Top)),] #Remove genes Gm- or -Rik or Ribosome and Mito genes
  Top5genes_cluster_D7[,i] = rownames(Top[1:5,])
}
colnames(Top5genes_cluster_D7) = c("Plasma.de.D7", "B.de.D7", "HNK.de.D7", "NK.de.D7", "Neut.de.D7", "Stromal.de.D7",
                                   "CD8T.de.D7", "CD4T.de.D7", "Il4i1DC.de.D7", "mac.de.D7", "Infmono.de.D7", "mono.de.D7",
                                   "CD103DC.de.D7", "CD11bDC.de.D7", "compmac.de.D7")


D9_markers_list = list(Plasma.de.D9, B.de.D9, HNK.de.D9, NK.de.D9, Neut.de.D9, Stromal.de.D9,
                       CD8T.de.D9, CD4T.de.D9, Il4i1DC.de.D9, mac.de.D9, Infmono.de.D9, mono.de.D9,
                       CD103DC.de.D9, CD11bDC.de.D9, compmac.de.D9)

Top5genes_cluster_D9 <- data.frame(matrix(ncol = 15, nrow = 5))
for (i in 1:length(D9_markers_list)){
  Top = D9_markers_list[[i]]
  Top$avg_log2FC = abs(Top$avg_log2FC)
  Top = Top[which(Top$p_val_adj < 0.05 | Top$avg_log2FC> 0.25),]
  Top = Top[ order(Top$p_val_adj,-Top$avg_log2FC),]
  Top = Top[!grepl("^Gm|Rik$|^Rps|^Rpl|^mt", rownames(Top)),] #Remove genes Gm- or -Rik or Ribosome and Mito genes
  Top5genes_cluster_D9[,i] = rownames(Top[1:5,])
}
colnames(Top5genes_cluster_D9) = c("Plasma.de.D9", "B.de.D9", "HNK.de.D9", "NK.de.D9", "Neut.de.D9", "Stromal.de.D9",
                                   "CD8T.de.D9", "CD4T.de.D9", "Il4i1DC.de.D9", "mac.de.D9", "Infmono.de.D9", "mono.de.D9",
                                   "CD103DC.de.D9", "CD11bDC.de.D9", "compmac.de.D9")

write.csv(Top5genes_cluster_D9,paste(path,"/Top5genes_cluster_D9.csv", sep = ""), row.names = T)
write.csv(Top5genes_cluster_D7,paste(path,"/Top5genes_cluster_D7.csv", sep = ""), row.names = T)

unique(Idents(TNBC))

# Set cell identity classes using SetIdent
cells.use <- WhichCells(TNBC, expression=`Cd8a` > 0 |`Cd8b1` > 0 , idents = c("8"))
TNBC <- SetIdent(TNBC, cells = cells.use, value = '2')
TNBC$seurat_clusters = Idents(TNBC)

TNBC <- RenameIdents(TNBC, '9' = '8')
TNBC$seurat_clusters = Idents(TNBC)

cells.use <- WhichCells(TNBC, expression=`Slamf7` > 0 , idents = c("Plasma Cells"))
TNBC <- RenameIdents(TNBC, '15' = '8')
TNBC$seurat_clusters = Idents(TNBC)

TNBC[[]]
Idents(TNBC) = "new.cluster.ids"
levels(Idents(TNBC))
library(MAST)
Plasma_Bcelldif =  FindMarkers(TNBC, ident.1 = "Plasmablast", ident.2 = "B Cell",
                               verbose = FALSE, test.use = "MAST")

CompMac_macdif =  FindMarkers(TNBC, ident.1 = "Complement Macrophage", 
                              ident.2 = "Macrophage",
                               verbose = FALSE, test.use = "MAST")

ChemMon_macdif =  FindMarkers(TNBC, ident.1 = "Chemoattractant Monocyte", 
                              ident.2 = "Monocytes",
                              verbose = FALSE, test.use = "MAST")

NKdif =  FindMarkers(TNBC, ident.1 = "Helper NK", 
                              ident.2 = "NK Cell",
                              verbose = FALSE, test.use = "MAST")

NKdif$gene = rownames(NKdif)
DoHeatmap(TNBC, 
          features = c(
            "Ptprc", #CD45
            "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Foxp3", # Tcells
            "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
            #"Mcpt1", "Fcer1a", "Kit", #Mast Cells
            "Stmn1", "Rrm1", #Helper NK/proliferating
            "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
            "Col4a1","Flt1","Tie1", "Col3a1", "Col1a1", "Dcn", #Stromal
            # Cd11b. cd11c.  dc sign.  cd103             f4/80
            "Itgam", "Itgax", "Cd209a", "Itgae","Adgre4","Adgre1","Cd80", "Clec9a", 
            "Il4i1", "Flt3", "Ly75",
            "Cd14", "Lyz1", #mon
            "Hp", "Mapk1",#chemoattractant
            "Mki67", "Top2a", # proliferation
            "C1qa", "C1qb","C1qc","C3", "C1ra", "C2", #complement
            "Cx3cr1", "Cd68","Ms4a7",#Mac
            
            "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2", "Il1a", "Il1b","Ccr2", #Mac
            "S100a8", "S100a9", "Ly6g","Mmp8", "Mmp9", #Neut
            "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
            
            "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
            "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Slamf7", "Cd44", "Prdm1", "Xbp1", "Mcm5"  #Plasma 
          ), slot = "scale.data", assay = "integrated",
          size = 3) + scale_fill_gradientn(colors = c("blue","grey","red")) + NoLegend()



DotPlot(TNBC, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Foxp3", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  #"Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Stmn1", "Rrm1", #Helper NK/proliferating
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Tie1", "Col3a1", "Col1a1", "Dcn", #Stromal
  # Cd11b. cd11c.  dc sign.  cd103             f4/80
  "Itgam", "Itgax", "Cd209a", "Itgae","Adgre4","Adgre1","Cd80", "Clec9a", 
  "Il4i1", "Flt3", "Ly75",
  "Cd14", "Lyz1", #mon
  "Hp", "Mapk1",#chemoattractant
  "Mki67", "Top2a", # proliferation
  "C1qa", "C1qb","C1qc","C3", "C1ra", "C2", #complement
  "Cx3cr1", "Cd68","Ms4a7",#Mac
  
  "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2", "Il1a", "Il1b","Ccr2", #Mac
  "S100a8", "S100a9", "Ly6g","Mmp8", "Mmp9", #Neut
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Slamf7", "Cd44", "Prdm1", "Xbp1", "Mcm5"  #Plasma 
)
, group.by = "new.cluster.ids", scale = T, assay = "RNA",
col.max = 8, col.min = 0
) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism(base_family = "Arial", base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1))



DotPlot(TNBC, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  #"Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  #"Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Nkx2-3", "Tie1", #Stromal
  # Cd11b. cd11c.  dc sign.  cd103             f4/80
  "Itgam", "Itgax", "Cd209a", "Itgae","Adgre4","Adgre1","Cd80", "Clec9a",
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #Myeloid MHC-II
  "C1qa", "C1qb","C1qc", "C3", "C1ra", "C2", #complement
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Ighe", "Slamf7", "Cd44", "Prdm1", "Xbp1", "Mcm5"  #Plasma 
)
, group.by = "new.cluster.ids", scale = F) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism(base_family = "Arial", base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1))

DotPlot(TNBC, group.by = "new.cluster.ids", features = c( "Ptprc", "Cd33",
                                           "Cd3e", "Cd4", "Cd8a","Ctla4", 
                                           "Ncr1", "Klrb1c", "Klre1", 
                                           "Cd79a", "Ms4a1", "Ighm", 
                                           "Cd300a", "Ly6c2", "Cd14", "Ccr2", 
                                           "Smox", "Apoe", "Ms4a7", "Ccr5", "Fcgr1", "Adgre1", "Cd68",
                                           "Ftl1", "Fth1", 
                                           "Tmem119", "Trem2", "Aif1", "Csf1r",
                                           "S100a9", "Camp", "Retnlg", "Ly6g",
                                           "Ifitm1", "Siglech", "Klk1", 
                                           "Zbtb46", "Itgax", "Cd86", "Cd209a",
                                           "Bst2", "Cmah", "Ly6a", "Nrp1", "Clec4g", 
                                           "Cacnb3", "Fscn1", "Syn3", "Tmem150c", 
                                           "Snca", "Ube2o", "Cd207","Clec9a", "Cd1d1", "Prdm1"
), scale = T) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

DotPlot(TNBC, group.by = "new.cluster.ids", features = c("Cd209a", "Btla", "Flt3", "Itgax", #cd11c
  "Clec9a", "Ly75", "Sirpa", "Bst2", "Itgae", #CD103
  "Itgam", #Cd11b, 
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa",
  "Cd80", "Cd86", "Il4i1"
), assay = "RNA", scale = F) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

plots <- VlnPlot(TNBC, features = c("Cd19", "Ms4a1", "Sdc1", "Tnfrsf13b", "Slamf7", "Prdm1", "Cd74", "Mcm5"), 
                 split.by = "sample",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)


plots <- VlnPlot(TNBC, features = c("Itgam", "Clec9a", "Cd209a", "Itgae", "Tnfrsf13b"), 
                 split.by = "sample", assay = "RNA",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

plots <- VlnPlot(TNBC, features = c("Ptprc", "Flt1", "Ly6a",
                                                   "Id3", "Fhl2", "Gng11", "Cd34", "Acvrl1",
                                                   "Itgb3", "Vcam1", "Thbd", "Cd55", "Vegfa", "Vegfd",
                                                   "Ece1", "Plec", "Ecscr","Tgfbr2"), split.by = "Antigen",
                 group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

plots <- VlnPlot(TNBC, features = c("Cd3e", "Cd4", "Cd8a","Cd8b1", "Ncr1", "Klre1", "Ptprc", "Klk1", "Nkg7", "Klrc2"), 
                 split.by = "sample",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)


plots <- VlnPlot(TNBC, features = c("Cd3e", "Cd4", "Cd8a","Cd8b1","Ighm"), 
                 split.by = "sample",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)
dev.off()

FeaturePlot(TNBC, 
            reduction = "umap", 
            features = c("Cd4", "Cd8a", "Cd8b1"), 
            #sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

Idents(TNBC) = "new.cluster.ids"
DimPlot(TNBC, label = T)
test <- WhichCells(TNBC, expression=`Cd8a` > 0 |`Cd8b1` >0 , idents = c("CD4 T"))
length(test)
test1 = WhichCells(TNBC, expression=`Cd4` > 0 , idents = c("CD4 T"))
length(test1)







DGE_allCells = FindAllMarkers(TNBC, log2FC.threshold = 0.25, test.use = "wilcox",
               min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
               assay = "RNA")

# DE for bulk deconvolution ----
TNBC[["CellTypes"]] = Idents(TNBC)
clusters.type = list(C1 = 'Plasma Cells', 
                     C2 = 'Fibroblast', 
                     C3 = c('Neutrophil', 'Mature B Cell', 'Immature DC'), 
                     C4 = c('NK Cell', 'Helper NK', 'TH17', 'CD4 T', 'CD8 T'),
                     C5 = c("Antigen-presenting B Cell","Chemoattractant Monocyte","Dendritic Cell",
                            "Macrophage","Monocytes","Complement Macrophage"))

cl.type = as.character(Idents(TNBC))

for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
unique(cl.type)
TNBC[["DeconvoCluster"]] = factor(cl.type, levels = names(clusters.type)) 
TNBC$DeconvoCluster
Idents(TNBC) = TNBC[["DeconvoCluster"]]
DeconvoMarkers <- FindAllMarkers(TNBC, log2FC.threshold = 0.2, test.use = "wilcox",
                                    min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                    assay = "RNA")
#write.csv(DeconvoMarkers,"/Users/lailarad/Documents/UMich Research/bulk deconvolution/DeconvoMarkers.csv",row.names = T)
#write.csv(DGE_allCells,"/Users/lailarad/Documents/UMich Research/bulk deconvolution/DeconvoAllMarkers.csv",row.names = T)


# pseudobulk ----
pseudobulk = AggregateExpression(
  TNBC,
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "new.cluster.ids",
  add.ident = NULL,
  slot = "data",
  verbose = TRUE
)

Idents(TNBC) = "Day"
pseudobulkD9 = AggregateExpression(
  subset(TNBC,idents = c("D9")),
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "new.cluster.ids",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)
Idents(TNBC) = "new.cluster.ids"
head(pseudobulkD9)

Idents(TNBC) = "Day"
pseudobulkD7 = AggregateExpression(
  subset(TNBC,idents = c("D7")),
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "new.cluster.ids",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)
Idents(TNBC) = "new.cluster.ids"
head(pseudobulkD7)


signature_genes1 = read.csv("/Users/lailarad/Documents/BI_EAE/aaron_pEAE/Signature_genes1.csv",row.names = NULL)
signature_genesTx = read.csv("/Users/lailarad/Documents/BI_EAE/Cohort1_RNA_seq/signature_genes_VLA4tx.csv",row.names = NULL)
signature_genes_day9_pEAE_bulk = read.csv("/Users/lailarad/Documents/BI_EAE/aaron_pEAE/aaronpEAE_signature_genes_day9.csv",row.names = NULL)


pseudobulkset = pseudobulkD9
mat = pseudobulkset$RNA
mat=log(mat + 1, base = 2)
mat = mat - rowMeans(mat)



sigGenes = (intersect(rownames(mat), signature_genes1$x))
sigGenes = (intersect(rownames(mat), signature_genes_day9_pEAE_bulk$x))
sigGenes = (intersect(rownames(mat), signature_genesTx$genes))

pheatmap::pheatmap(t(mat[sigGenes,]),
                   #annotation_col = anno,
                   show_colnames = T,
                   fontsize_row = 18,
                   fontsize_col = 12,
                   cluster_cols = T,
                   cluster_rows = T,
                   #annotation_colors = ann_colors,
                   #filename = paste(path, "signatureCelltypesD79_vla4tx_log2.png",sep=""),
                   width = 8,
                   height = 5
                   
                   
)
dev.off()
getwd()
plots <- VlnPlot(TNBC, features = sigGenes[(25:31)], 
                 split.by = "sample",
                 group.by = "new.cluster.ids", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 3)

DotPlot(TNBC, features = sigGenes) + 
  RotatedAxis()

#pseudobulk conditions = sample----
TNBC[[]]
Idents(TNBC) = "sample"
pseudobulkSamples = AggregateExpression(
  TNBC,
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "sample",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)

head(pseudobulkSamples)

library(ggbiplot)
library(RColorBrewer)
library(caret)
plot_biplot = function(data, genesEN){
  data = pseudobulkSamples$RNA
  genesEN = sigGenes
  g_EN_value <- unique(rownames(data)[unlist(genesEN)])
  n_ENvalue <- na.omit(data[sigGenes, ])
  ggbiplot(prcomp(t(n_ENvalue), scale.=T), ellipse=T, var.axes=F, 
             var.scale=1, circle=T) + theme_classic() + 
      #ggtitle(paste("EAE Anti-VLA4 alpha=", names(genesEN$ENgenes[num]),"Geneset")) + 
      geom_point(size=3) +
      scale_color_manual(name="Group")
}


    ggbiplot(prcomp(t(pseudobulkSamples$RNA[sigGenes,]), scale.=T), ellipse=T,  var.axes=F, 
             var.scale=1, circle=T) + theme_classic() + 
      #ggtitle(paste("EAE Anti-VLA4 alpha=", names(genesEN$ENgenes[num]),"Geneset")) + 
      geom_point(size=3, color=colSide) +
      scale_color_manual(name="Group", values = ecolor)
  

sigGenes = (intersect(rownames(pseudobulkSamples$RNA), signature_genes1$x))
plot_biplot(pseudobulkSamples, sigGenes)
ggbiplot(prcomp(t(sigGenes), scale.=T), ellipse=T, groups=names(colSide), var.axes=F, 
         var.scale=1, circle=T) + theme_classic() + 
  ggtitle(paste("EAE Anti-VLA4 alpha=", names(genesEN$ENgenes[num]),"Geneset")) + 
  geom_point(size=3, color=colSide) +
  scale_color_manual(name="Group", values = ecolor)

sample_info3 = sample_info[-which(sample_info$Immunization == "OVA"),]
groups_t4 = factor(colnames(pseudobulkSamples$RNA))
coul <- colorRampPalette(brewer.pal(11, "PuOr"))(4)
colSide = as.character(colnames(pseudobulkSamples$RNA))
Conditions = sort(unique(colnames(pseudobulkSamples$RNA)))

for (x in 1:length(groups_t4)) {
  if (groups_t4[x] == Conditions[1]) {colSide[x] <- "black" #purple = healthy
  } else if (groups_t4[x] == Conditions[2]) {colSide[x] <- "#6a329f" #  yellow = "MOG _ Anti-VLA4 _ Day 10" 
  } else if (groups_t4[x] == Conditions[3]) {colSide[x] <- "#ff7f00" #orange = "MOG _ Control Ab _ Day 10"
  } else if (groups_t4[x] == Conditions[4]) {colSide[x] <- "#e31a1c" #red = "MOG _ Pre-Tx"
  } }


names(colSide) <- groups_t4
ecolor = c("black", 
           "#6a329f",
           "#ff7f00",
           "#e31a1c"
)


count.data = pseudobulkSamples$RNA
is.na(count.data) %>% table() # Check for NA values
count.data[is.na(count.data)] <- 0 # Convert NA to 0
is.na(count.data) %>% table() # Check again for NA values
count.data <- count.data[rowSums(count.data)>0,] # Remove 0 expression genes
sample_num = length(count.data)
count.data <- count.data[rowSums(count.data == 0) <= sample_num*(3/4),] # Remove mostly 0 genes
count.data <- count.data[rowSums(count.data) > sample_num, ] # Remove low expression genes
count.data = count.data[!grepl("^Gm|Rik$", rownames(count.data)),] #Remove genes Gm- or -Rik


xfactors_eae = cbind.data.frame(groups_t4)
rownames(xfactors_eae)
dataset_eae = t(count.data )
dataset_eae = dataset_eae[(rownames(dataset_eae) %in% groups_t4),]
dataset_eae = as.data.frame(cbind(xfactors_eae, dataset_eae))
dataset_eae <- dataset_eae[complete.cases(dataset_eae), ]
set.seed(123) 
training.samples <- dataset_eae$groups_t4 %>% createDataPartition(p = 1.0, list = F)
train.data  <- dataset_eae[training.samples, ]
test.data <- dataset_eae[-training.samples, ]



x <- model.matrix(groups_t4~., train.data)[,-1] # Predictor variables
y <- train.data$groups_t4 # Outcome variable

multinom4_EN <- function(x, y, num_folds){
  alphaList <- seq(1.0, 0.0, -0.05)
  pb <- txtProgressBar(min = 0, max = length(alphaList), style = 3)
  ENgenes <- lapply(1:length(alphaList), function(i) {
    setTxtProgressBar(pb, i)
    # cvfit <- suppressMessages(cv.glmnet(x, y, alpha = alphaList[i], 
    #                                     nfolds = num_folds, family = "multinomial"))
    # best.lambda <- cvfit$lambda.min
    fit <- (glmnet(x, y, alpha = alphaList[i], #lambda = best.lambda, 
                   family = "multinomial"))
    Coef<-coef(fit, s = best.lambda)
    Index<-c(Coef[["0"]]@i[-1], Coef[["1"]]@i[-1], Coef[["2"]]@i[-1], Coef[["3"]]@i[-1])
  })
  
  names(ENgenes) <- alphaList
  return(list("ENgenes"=ENgenes))
  
}

genes_EN_eae4 <- multinom4_EN( x, y, 100)

#bootstrap genes for each sample
#https://bootcamp.biostars.io/archives/2016/day3/docs/BootstrapRNAseq.html

samp=c(sample(1,6,replace=TRUE),sample(1,6,replace=TRUE), sample(1,6,replace=TRUE), sample(1,6,replace=TRUE))
ReadCounts = pseudobulkSamples$RNA
bReadCounts=ReadCounts[,samp]
colnames(bReadCounts)=c(rep(colnames(ReadCounts)[1],6),rep(colnames(ReadCounts)[2],6),
                        rep(colnames(ReadCounts)[3],6), rep(colnames(ReadCounts)[4],6))

ReadCountsAdj=as.matrix(ReadCounts)
ReadCountsAdj=ReadCountsAdj-.3
ReadCountsAdj[ReadCountsAdj<0]=0.25
nn=nrow(ReadCountsAdj)*ncol(ReadCountsAdj)
bNoisyReads=matrix(rpois(nn,ReadCountsAdj),ncol=ncol(ReadCountsAdj))
colSums(ReadCounts)
colSums(bNoisyReads)


ReadCounts0=as.matrix(ReadCounts)
ReadCounts0[ReadCounts0==0]=.25
libSizes=colSums(ReadCounts0)
ReadProp=ReadCounts0
for (i in 1:ncol(ReadCounts0)) ReadProp[,i]=ReadCounts0[,i]/libSizes[i]
means=rowMeans(ReadProp)
SDs=apply(ReadProp,1,sd)
oGenes=order(means,decreasing=TRUE)
median(means)
head(means[oGenes[1:10]])
median(SDs)
head(SDs[oGenes[1:10]])

bProp=ReadProp
bMeans=bProp
for (i in 1:nrow(ReadProp)) bProp[i,]=rnorm(ncol(ReadProp),means[i],SDs[i])
for (i in 1:ncol(bProp))  bMeans[,i]=bProp[,i]*libSizes[i]
bMeans[bMeans<=0.25]=0.25
bReads=matrix(rpois(nn,bMeans),ncol=ncol(bMeans))
colnames(bReads)=colnames(ReadCounts)
head(ReadCounts)
head(bReads)

# bootstrap single cells to make pseudobulk into multiple samples ----


#DEG for cell types ----
#Tcells ----
Idents(TNBC)
Tcell_data = subset(TNBC,idents = c("CD4 T", "CD8 T"))
CD4Tcell_data = subset(TNBC,idents = c("CD4 T"))
CD8Tcell_data = subset(TNBC,idents = c("CD8 T"))
Idents(CD4Tcell_data) = "sample"
CD4_D7_OvM <- FindMarkers(CD4Tcell_data, ident.1 = "D7 MOG", ident.2 = "D7 OVA")
CD4_D9_OvM <- FindMarkers(CD4Tcell_data, ident.1 = "D9 MOG", ident.2 = "D9 OVA")

Idents(CD8Tcell_data) = "sample"
CD8_D7_OvM <- FindMarkers(CD8Tcell_data, ident.1 = "D7 MOG", ident.2 = "D7 OVA")
CD8_D9_OvM <- FindMarkers(CD8Tcell_data, ident.1 = "D9 MOG", ident.2 = "D9 OVA")

library(AnnotationHub)
library(ensembldb)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

  
## feature 1: numeric vector
geneList = CD8_D7_OvM$avg_log2FC
geneList =  CD4_D9_OvM$avg_log2FC
geneList = CD4_D7_OvM$avg_log2FC


## feature 2: named vector
names(geneList) = as.character(rownames(CD8_D7_OvM))
names(geneList) = as.character(rownames(CD4_D9_OvM))
names(geneList) = as.character(rownames(CD4_D7_OvM))

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)

geneList <- geneList[!is.na(geneList)]

gse <- gseGO(geneList=geneList, 
             ont ="BP", 
             keyType = "SYMBOL", 
             #nPerm = 10000, 
             minGSSize = 10, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.1, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")
dotplot(gse, showCategory=50, split=".sign",
        x = "GeneRatio",
        label_format = function(x) stringr::str_wrap(x, width=100),
        title = "Mog9 vs. OVA") + 
  facet_grid(.~.sign) +
  #theme_prism(base_family = "Arial", base_size = 16) +
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  ) 
gse1 = gse
gse.df = as.data.frame(gse)
filter(gse, geneList == "Cyp2e1" )
gse@result$core_enrichment = list

class(gse1@result$core_enrichment)

heatplot(gse, foldChange=geneList, showCategory=50,# symbol = "dot",
         label_format = function(x) stringr::str_wrap(x, width=50)) +
  #theme_prism(base_family = "Arial", base_size = 16) +
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  ) 
  

#new T cell subclusters ----
getwd()
TsubClusters = readRDS("newTcellclusters.rds")
Idents(TsubClusters)
TsubClusters = RenameIdents(TsubClusters, c("Ferritin CD3+ Cells" = "Ferritin Myeloid Cells",
                                            "Cytotoxic CD8" = "Cytotoxic CD8 Tem",
                                            "CD8 Tcm" = "Naive CD8"
                                            ))
TsubClusters$CellTypes = Idents(TsubClusters) 
table(TsubClusters$CellTypes)

DimPlot(TsubClusters, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T)

getPalette = brewer.pal(n = 8, name = "Dark2")

png(file = "UMAPnewsubclusterTcells.png",
    units = "in",width = 11, height = 11, res = 400)
DimPlot(TsubClusters, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T, cols = getPalette) + 
  theme_prism(base_family = "Arial", base_size = 24) + theme(legend.text = element_text(size = 28)) +
  NoLegend()
dev.off()

Tcellmarkers = c(
  "Cd3e", "Cd8a", "Cd4",
  #naive T cells:
  "Ccr7", "Cd28",
  #cytotoxic CD8:
  "Prf1","Gzmb","Gzma","Ccl4","Ccl5",
  #Th1:
  "Ccr4","Cxcr3", "Ccr5", "Il12rb1","Il12rb2", "Ifng","Ifngr1", "Tbx21", #Tbet 
  #Th2:
  "Il4","Il4ra", "Il5","Il5ra", "Gata3",
  #Th9:
  "Ccr3", "Ccr6", "Spi1","Il9r",#"Il22ra1",
  #Th17: Ccr6, Ccr4 Nk1.1
  "Klrb1c", "Il17a","Il17f","Il6ra","Tgfb1","Il23a",
  #Treg: CD127 = Il7r
  "Il7r", "Il2ra", "Ctla4", "Il2", "Foxp3", "Il10","Il10ra",
  #Tfh:
  "Cxcr5", "Cd40lg", "Icos", "Bcl6", "Il21r",
  #gammadelta T cell:
  "Il23r", "Rorc","Rora","Tcrg-V6",
  #Memory:
  "Cd44","Sell", #CD62L
  "Cd27", "Itgae", #Cd103
  "Lag3",
  "Izumo1r", #Fr-4
  "Nt5e", #CD73
  "Slamf6", "Xcl1","Cx3cr1", "Tox", "Gzmk", "Pdcd1",
  "Trgv2", "Trdc", "Trdv4",
  "Klk1", "Nkg7", "Klrc2", "Ncr1", "Klre1", "Itgam","Itgax",
  "Ftl1", "Fth1", "Apoe"
)
getwd()
png(file = "Tcell_clusterbyhanddotplot.png",
    units = "in",width = 20, height = 10, res = 400)
col_grid <- rgb(235, 235, 235, 100, maxColorValue = 255)
DotPlot(TsubClusters, features = Tcellmarkers,
        scale = F#, cols = c("blue", "red")
) + RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +
  theme(panel.grid.major.x = element_line(colour = col_grid),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
dev.off()

BiocManager::install("pals")
library("pals")

# Frequency plots
TCellTypeTable = as.data.frame(proportions(table(Idents(TsubClusters),TsubClusters@meta.data$sample), margin = 2))
TCellTypeTable$Freq = TCellTypeTable$Freq*100
levels(TCellTypeTable$Var1)


#DiscretePalette(7, palette = "watlington", shuffle = F)
#as.vector(cols25(7))
getPalette
# Discrete
#pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, okabe, polychrome, stepped, stepped2, stepped3, tol, watlington,
#          main="Discrete", show.names=FALSE)
#pal.bands(cols25)


png(file = "TCellProportionsBMES.png",
    units = "in",width = 7, height = 5, res = 400)
ggplot(TCellTypeTable,aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism(base_family = "Arial", base_size = 16)+
  labs(x = "Condition", y= "Percent of Cells") + 
  scale_fill_manual(values = getPalette, name = "Cell Type") + 
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  )
dev.off()

TCellTypeData = TCellTypeTable
# load stringr library
library(stringr)

# Split name column into firstname and last name
TCellTypeData[c('Day', 'Antigen')] <- str_split_fixed(TCellTypeData$Var2, ' ', 2)

# plotCellType = function(data, celltype){
#   ggplot(data[which(data$Var1 == celltype),], aes(x = Antigen, y= Freq, fill = Day)) + 
#     geom_bar(position = position_dodge(),stat = 'identity') +
#     scale_y_continuous(expand = expansion(mult = c(0, .1)))+
#     theme_prism(base_family = "Arial", base_size = 5)+
#     labs(x = "Condition", y= "Percent of Cells", title = paste(celltype)) + 
#     theme(legend.text = element_text(size = 5)
#     )
# }

myplots = lapply(unique(TCellTypeData$Var1), FUN = plotCellType, data = TCellTypeData)

ggarrange(plotlist = myplots, ncol=4, nrow=2, 
          common.legend = T, legend = "right")
dev.off()

plotCellTypePresentation = function(data,celltype){
  ggplot(data[which(data$Var1 == celltype),], aes(x = Antigen, y= Freq, fill = Day)) + 
    geom_bar(position = position_dodge(),stat = 'identity') +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    theme_prism(base_family = "Arial", base_size = 16)+
    labs(y= "Percent of Cells", title = paste(celltype)) + 
    theme(legend.text = element_text(size = 16),
          axis.title.x = element_blank()
    )
}

plotCellTypePresentationCounts = function(data,celltype){
  ggplot(data[which(data$Var1 == celltype),], aes(x = Antigen, y= Freq, fill = Day)) + 
    geom_bar(position = position_dodge(),stat = 'identity') +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))+
    theme_prism(base_family = "Arial", base_size = 16)+
    labs(y= "Cell Counts", title = paste(celltype)) + 
    theme(legend.text = element_text(size = 16),
          axis.title.x = element_blank() 
    )
}

png(file = "NaiveCD8_per_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentation(TCellTypeData, "Naive CD8")
dev.off()
png(file = "NaiveCD8_count_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentationCounts(TCellTypeCounts, "Naive CD8")
dev.off()

png(file = "TemCD8_per_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentation(TCellTypeData, "Cytotoxic CD8 Tem")
dev.off()
png(file = "TemCD8_count_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentationCounts(TCellTypeCounts, "Cytotoxic CD8 Tem")
dev.off()

png(file = "NaiveCD4_per_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentation(TCellTypeData, "Naive CD4")
dev.off()
png(file = "NaiveCD4_count_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentationCounts(TCellTypeCounts, "Naive CD4")
dev.off()

png(file = "TemCD4_per_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentation(TCellTypeData, "CD4+ Effector Memory")
dev.off()
png(file = "TemCD4_count_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentationCounts(TCellTypeCounts, "CD4+ Effector Memory")
dev.off()

getwd()
png(file = "Treg_per_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentation(TCellTypeData, "Treg")
dev.off()

png(file = "Treg_count_BMES.png",units = "in",width = 5, height = 5, res = 400)
plotCellTypePresentationCounts(TCellTypeCounts, "Treg")
dev.off()


plotCellTypePresentation(CellTypeData, "CD8 T")
plotCellTypePresentation(CellTypeData, "Plasmablast")

myplots = lapply(c("B Cell", "CD8 T", "Plasmablast"), FUN = plotCellTypePresentation, data = CellTypeData)

#png(file = "ImmuneCell_prop_BMES.png",units = "in",width = 16, height = 5, res = 400)
ggarrange(plotlist = myplots, ncol=3, nrow=1, 
          common.legend = T, legend = "right")
dev.off()


TCellTypeCounts =as.data.frame(table(Idents(TsubClusters),TsubClusters@meta.data$sample))

# Split name column into firstname and last name
TCellTypeCounts[c('Day', 'Antigen')] <- str_split_fixed(TCellTypeCounts$Var2, ' ', 2)

myplots = lapply(unique(TCellTypeData$Var1), FUN = plotCellType, data = TCellTypeCounts)

ggarrange(plotlist = myplots, ncol=4, nrow=4, 
          common.legend = T, legend = "right")
dev.off()

CD8Pathogenicmarkers = c(
"CD86", "Il2rb", #CD122-CD8+ T cells pathogenic
"Il17f",
"Rorc", 
"Ifng", "Klrb1c", "Tnf", 
"Cd7", "Ccl5", "Csf2",
)

TemCD8_D9_OvM$gene[1:7]

plots <- VlnPlot(TsubClusters, features = c("Il1b","Il23a","Il23r", "Il1r1","Itga4","Itgal",
                                            "Mmp2", "Mmp8", "Mmp9"
                                            ), 
group.by = "CellTypes", 
split.by = "sample",pt.size = 0, combine = FALSE) 
wrap_plots(plots = plots, ncol = 3)

png(file = "CD8Tcell_PathogenicMarkers.png",units = "in",width = 16, height = 15, res = 400)
wrap_plots(plots = plots, ncol = 3)
dev.off()


#DEG Functional Enrichment new subcluster ----
#Tcells ----
Idents(TsubClusters)

NaiveCD4 = subset(TsubClusters,idents = c("Naive CD4"))
NaiveCD8 = subset(TsubClusters,idents = c("Naive CD8"))
TemCD4 = subset(TsubClusters,idents = c("CD4+ Effector Memory"))
TemCD8 = subset(TsubClusters,idents = c("Cytotoxic CD8 Tem"))
Treg = subset(TsubClusters,idents = c("Treg"))

Idents(Treg)
head(Treg[[]])
pseudobulk_Treg = AggregateExpression(
  Treg,
  assays = "RNA",
  features = NULL,
  return.seurat = FALSE,
  group.by = "sample",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)
head(pseudobulk_Treg)
class(pseudobulk_Treg$RNA)
write.csv(pseudobulk_Treg$RNA,paste0(path,"EAE_scRNAseq_Treg_pseudobulkbysample.csv"), row.names = T)


#Find Markers: Positive values indicate that the gene is more 
#highly expressed in the first group


Idents(NaiveCD4) = "sample"
NCD4_D7_OvM <- FindMarkers(NaiveCD4, ident.1 = "D7 MOG", ident.2 = "D7 OVA")
NCD4_D9_OvM <- FindMarkers(NaiveCD4, ident.1 = "D9 MOG", ident.2 = "D9 OVA")

Idents(NaiveCD8) = "sample"
NCD8_D7_OvM <- FindMarkers(NaiveCD8, ident.1 = "D7 MOG", ident.2 = "D7 OVA")
NCD8_D9_OvM <- FindMarkers(NaiveCD8, ident.1 = "D7 MOG", ident.2 = "D7 OVA")

Idents(TemCD4) = "sample"
TemCD4_D7_OvM <- FindMarkers(TemCD4,ident.1 = "D7 OVA", ident.2 = "D7 MOG")
TemCD4_D9_OvM <- FindMarkers(TemCD4, ident.1 =  "D7 OVA", ident.2 = "D7 MOG")

Idents(TemCD8) = "sample"
TemCD8_D7_OvM <- FindMarkers(TemCD8, ident.1 = "D7 MOG", ident.2 = "D7 OVA")
TemCD8_D9_OvM <- FindMarkers(TemCD8, ident.1 = "D9 MOG", ident.2 = "D9 OVA")

Idents(Treg) = "sample"
Treg_D7_OvM <- FindMarkers(NaiveCD4, ident.1 = "D7 OVA", ident.2 = "D7 MOG", min.pct = 0.001)
Treg_D9_OvM <- FindMarkers(NaiveCD4, ident.1 = "D7 OVA", ident.2 = "D7 MOG", min.pct = 0.01)


library(AnnotationHub)
library(ensembldb)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)



TemCD4_D9_OvM$gene = rownames(TemCD4_D9_OvM)
TemCD4_D9_OvM = inner_join(TemCD4_D9_OvM, annotations, by=c("gene"="gene_name")) 
TemCD4_D7_OvM$gene = rownames(TemCD4_D7_OvM)
TemCD4_D7_OvM = inner_join(TemCD4_D7_OvM, annotations, by=c("gene"="gene_name")) 

TemCD8_D9_OvM$gene = rownames(TemCD8_D9_OvM)
TemCD8_D9_OvM = inner_join(TemCD8_D9_OvM, annotations, by=c("gene"="gene_name")) 
TemCD8_D7_OvM$gene = rownames(TemCD8_D7_OvM)
TemCD8_D7_OvM = inner_join(TemCD8_D7_OvM, annotations, by=c("gene"="gene_name")) 


Treg_D9_OvM$gene = rownames(Treg_D9_OvM)
Treg_D9_OvM = inner_join(Treg_D9_OvM, annotations, by=c("gene"="gene_name")) 
Treg_D7_OvM$gene = rownames(Treg_D7_OvM)
Treg_D7_OvM = inner_join(Treg_D7_OvM, annotations, by=c("gene"="gene_name")) 



## feature 1: numeric vector
geneList = NCD8_D7_OvM$avg_log2FC
geneList = NCD8_D9_OvM$avg_log2FC
geneList =  NCD4_D9_OvM$avg_log2FC
geneList = NCD4_D7_OvM$avg_log2FC

geneList1 = TemCD8_D7_OvM$avg_log2FC
geneList2 = TemCD8_D9_OvM$avg_log2FC
geneList4 =  TemCD4_D9_OvM$avg_log2FC
geneList3 = TemCD4_D7_OvM$avg_log2FC

geneList5 = Treg_D9_OvM$avg_log2FC


## feature 2: named vector
names(geneList) = as.character(rownames(NCD8_D7_OvM))
names(geneList) = as.character(rownames(NCD8_D9_OvM))
names(geneList) = as.character(rownames(NCD4_D9_OvM))
names(geneList) = as.character(rownames(NCD4_D7_OvM)) #nothing enriched

names(geneList1) = as.character((TemCD8_D7_OvM$entrezid)) #nothing enriched
names(geneList2) = as.character((TemCD8_D9_OvM$entrezid))
names(geneList4) = as.character(TemCD4_D9_OvM$entrezid)
names(geneList3) = as.character(TemCD4_D7_OvM$entrezid)
names(geneList5) = as.character(Treg_D7_OvM$entrezid)

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)
geneList1 = sort(geneList1, decreasing = TRUE)
geneList2 = sort(geneList2, decreasing = TRUE)
geneList3 = sort(geneList3, decreasing = TRUE)
geneList4 = sort(geneList4, decreasing = TRUE)
geneList5 = sort(geneList5, decreasing = TRUE)

geneList <- geneList[!is.na(geneList)]
geneList1 <- geneList1[!is.na(geneList1)]
geneList2 <- geneList2[!is.na(geneList2)]
geneList3 <- geneList3[!is.na(geneList3)]
geneList4 <- geneList4[!is.na(geneList4)]
geneList5 <- names(geneList5)[!is.na(names(geneList5))]
geneList5 <- geneList5[!grepl("^NA|^c", names(geneList5))]
CD4MOG7v9[!grepl("^Rp|^Rps|^H", rownames(CD4MOG7v9)),]

# install.packages("devtools")
devtools::install_github("immunogenomics/presto")
library(presto)

data(exprs)
head(exprs)[, 1:10]
data(y)
head(y)
length(y)
length(exprs[1,])

class(y)
class(exp_1)
exp_1 = matrix(c(2,4,1, 2, 1,2, 3, 4), nrow = 2, ncol = 4,
                dimnames = list(c("row1", "row2"),
                                c("C.1", "C.2", "C.3", "C.4")))
exp_2 = as.data.frame(exp_1)

factors = as.factor(c("A", "A", "B", "B"))
auc.res = wilcoxauc(exp_2, factors, groups_use = c('A', 'B'))

x = Treg@assays$RNA@counts
Idents(Treg) = "sample"
FM.Treg.genes.res = FindMarkers(object = Treg, ident.1 = "D9 OVA", ident.2 = "D9 MOG", 
                                min.pct = -Inf, logfc.threshold = -Inf, 
                                min.cells.feature = 1, min.cells.group = 1)

FM.Treg.genes.res$feature = rownames(FM.Treg.genes.res)

Treg.genes.res = (wilcoxauc(Treg, "sample", seurat_assay = 'RNA', groups_use = c('D7 OVA', 'D7 MOG')))
dplyr::count(Treg.genes.res, group)

FM.Treg.genes.res = inner_join(FM.Treg.genes.res, annotations, by=c("feature"="gene_name"))

Treg_geneList =  Treg.genes.res %>%
  arrange(desc(logFC)) %>%
  dplyr::select(entrezid, logFC)

Treg_geneList =  FM.Treg.genes.res %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(entrezid, avg_log2FC)

Tregccgenes$sign = "Up"
Tregccgenes[which(Tregccgenes$logFc < 0),]$sign = "Down"

head(Treg.genes.res1)
Treg.genes.res1 =  FM.Treg.genes.res
Treg.genes.res1$Day = "Day 7"
Treg.genes.res1$sign = "Up"
Treg.genes.res1[which(Treg.genes.res1$avg_log2FC < 0),]$sign = "Down"

  
  xx <- clusterProfiler::compareCluster(entrezid~Day+sign, data=Treg.genes.res1, fun = enricher,
                                        TERM2GENE=wp[,c("wpid", "gene")], TERM2NAME=wp[,c("wpid", "name")], 
                                        pvalueCutoff = 0.2)
getwd()
#png(file = "Treg_WPenrich_updowninMOG_FM.png",units = "in",width = 6, height = 4, res = 400)
clusterProfiler::dotplot(xx, x="Day",showCategory = 29,
                         title = "Treg Cells\nWikiPathways Enrichment Analysis",
                         label_format = function(x) stringr::str_wrap(x, width=60)) + 
  facet_grid(~sign) +
  #aes(x=fct_relevel(time, c('0h', '2h', '6h', '24h'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y =element_text(size = 12, face = "bold")) 

dev.off()



library(tibble)


convertMouseMGItoHumanSym <- function(x){
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org/")
  
  host="https://dec2021.archive.ensembl.org/"
  
  genesV2 <- getLDS(attributesL = c("hgnc_symbol"),
                    filters = "mgi_symbol",
                    values = x,
                    mart = mouse,
                    attributes = c("mgi_symbol"),
                    martL = human,
                    uniqueRows=TRUE)
  
  # resulting table is in a different order to the input list
  # reorder to get the output the right way around
  human_genes <- genesV2
  return(human_genes)
}
hsym = convertMouseMGItoHumanSym(FM.Treg.genes.res$symbol)
head(hsym)

Treg_ranks<- deframe(Treg_geneList)



head(ranks)

TsubClusters[[]]
wilcoxauc(seurat_object, 'group_name')


geneList = geneList5
data(geneList, package = "DOSE")

x = enrichPathway (gene = names(geneList5),
                   organism = "mouse")
head(x)

library(ReactomePA)
library(DOSE)
library(forcats)
library(msigdbr)
y = gsePathway(Treg_ranks,
               pvalueCutoff = 0.2,
               pAdjustMethod = "BH",
               organism = "mouse"
               )

all_gene_sets = msigdbr(species = "Mus musculus")
unique(all_gene_sets$gs_subcat)

m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

CTD_MS_pathways = read.csv( "/Users/lailarad/Downloads/CTD_disease_pathways_inferred_1719717436116.csv", check.names=FALSE,header=T, sep=",")
view(CTD_MS_pathways)

# Query AnnotationHub
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))
human_ens <- human_ens[["AH113665"]]
annotations_h <- genes(human_ens, return.type = "data.frame")

CTD_MS_pathways = inner_join(CTD_MS_pathways, annotations_h, by=c("InferenceGeneSymbol"="gene_name"))

convertMouseEntrezIDtoGene <- function(x){
  require("biomaRt")
  mouse1 <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                    verbose = TRUE, host = "dec2021.archive.ensembl.org")
  mouse2 <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                    verbose = TRUE, host = "dec2021.archive.ensembl.org")
  
  genesV2 <- getLDS(attributes = c("entrezgene_id"),
                    filters = "entrezgene_id",
                    values = x,
                    mart = mouse1,
                    attributesL = c("mgi_symbol"),
                    martL = mouse2,
                    uniqueRows=TRUE)
  
  # resulting table is in a different order to the input list
  # reorder to get the output the right way around
  mID_genes <- genesV2
  return(mID_genes)
}



FM.Treg.genes.res = FindMarkers(object = Treg, ident.1 = "D7 OVA", ident.2 = "D7 MOG", 
                                min.pct = -Inf, logfc.threshold = -Inf, 
                                min.cells.feature = 1, min.cells.group = 1)

FM.Treg.genes.res$symbol = rownames(FM.Treg.genes.res)


FM.Treg.genes.res = inner_join(FM.Treg.genes.res, annotations, by=c("feature"="gene_name"))


h_id = convertMouseMGItoHumanEntrez(FM.Treg.genes.res$symbol)
head(h_id)

colnames(h_id)[1] = "symbol"

FM2 = merge(FM.Treg.genes.res, h_id, by = c("symbol"))

Treg_geneList =  FM2 %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(NCBI.gene..formerly.Entrezgene..ID, avg_log2FC)

Treg_ranks<- deframe(Treg_geneList)

biocarta_gene_sets = msigdbr(species = "human", category = "C2", subcategory = "CP:BIOCARTA")


CTD_MS_pathways= CTD_MS_pathways %>% dplyr::select(PathwayName, entrezid)
colnames(CTD_MS_pathways) = c("gs_name", "entrez_gene")
CTD_MS_pathways$entrez_gene = as.integer(CTD_MS_pathways$entrez_gene)

CTD_gene_sets = data.frame(gs_name = CTD_MS_pathways$gs_name,entrez_gene = CTD_MS_pathways$entrez_gene, 
                           row.names = NULL)
y1 = GSEA(Treg_ranks,
         pvalueCutoff = 0.7,
         pAdjustMethod = "BH",
         #OrgDb = org.Mm.eg.db, 
         #TERM2GENE = CTD_gene_sets,
         TERM2GENE = biocarta_gene_sets[,c("gs_description", "entrez_gene")],
         minGSSize = 1
         
)

length(intersect(as.character(CTD_gene_sets$entrez_gene), names(Treg_ranks)))

head(y1)

y_filt1 = y1[which(abs(y1$NES)>0.2),]
y_filt1 = y_filt1[which((y_filt1$qvalue)< 0.75),]

y_filt1$Condition = "Downregulated"
y_filt1[which((y_filt1$NES)>0),]$Condition = "Upregulated"
y_filt1$Day = "Day 7"

FM.Treg.genes.res = FindMarkers(object = Treg, ident.1 = "D9 OVA", ident.2 = "D9 MOG", 
                                min.pct = -Inf, logfc.threshold = -Inf, 
                                min.cells.feature = 1, min.cells.group = 1)

FM.Treg.genes.res$symbol = rownames(FM.Treg.genes.res)


FM.Treg.genes.res = inner_join(FM.Treg.genes.res, annotations, by=c("feature"="gene_name"))

convertMouseMGItoHumanEntrez <- function(x){
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org/")
  
  host="https://dec2021.archive.ensembl.org/"
  
  genesV2 <- getLDS(attributesL = c("entrezgene_id"),
                    filters = "mgi_symbol",
                    values = x,
                    mart = mouse,
                    attributes = c("mgi_symbol"),
                    martL = human,
                    uniqueRows=TRUE)
  
  # resulting table is in a different order to the input list
  # reorder to get the output the right way around
  human_genes <- genesV2
  return(human_genes)
}
h_id = convertMouseMGItoHumanEntrez(FM.Treg.genes.res$symbol)
head(h_id)

colnames(h_id)[1] = "symbol"

FM2 = merge(FM.Treg.genes.res, h_id, by = c("symbol"))

Treg_geneList =  FM2 %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::select(NCBI.gene..formerly.Entrezgene..ID, avg_log2FC)

Treg_ranks<- deframe(Treg_geneList)



CTD_MS_pathways= CTD_MS_pathways %>% dplyr::select(PathwayName, entrezid)
colnames(CTD_MS_pathways) = c("gs_name", "entrez_gene")
CTD_MS_pathways$entrez_gene = as.integer(CTD_MS_pathways$entrez_gene)

CTD_gene_sets = data.frame(gs_name = CTD_MS_pathways$gs_name,entrez_gene = CTD_MS_pathways$entrez_gene, 
                              row.names = NULL)
y = GSEA(Treg_ranks,
         pvalueCutoff = 0.7,
         pAdjustMethod = "BH",
         #OrgDb = org.Mm.eg.db, 
         #TERM2GENE = CTD_gene_sets,
         TERM2GENE = biocarta_gene_sets[,c("gs_description", "entrez_gene")],
         minGSSize = 1
                
)

length(intersect(as.character(CTD_gene_sets$entrez_gene), names(Treg_ranks)))

head(y)

y_filt = y[which(abs(y$NES)>0.5),]
y_filt = y_filt[which((y_filt$qvalue)< 0.75),]

y_filt$Condition = "Downregulated"
y_filt[which((y_filt$NES)>0),]$Condition = "Upregulated"
y_filt$Day = "Day 9"

y_filt2 = rbind(y_filt, y_filt1)
y_filt2 = y_filt2[which((y_filt2$qvalue)<0.55),]

length(unique(CTD_gene_sets$gs_name))

png(file = "Treg_CTD_MSpathways_GSEA.png",units = "in",width = 18, height = 15, res = 400)
ggplot(y_filt2, aes(x = Day, y = fct_reorder(Description, NES))) + 
  geom_point(aes(color = NES, size = setSize) ) +
  theme_bw(base_size = 14) +
  scale_colour_gradient2( low = "blue",high="red", mid = "white",na.value = "black") +
  ylab(NULL) +
  ggtitle("Multiple Sclerosis Pathways GSEA from Tregs at the Scaffold") +
  facet_grid(~Condition) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 16),
        plot.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y =element_text(size = 16, face = "bold"),
        axis.title.x = element_blank(),
        strip.text.x  = element_text(size = 16, face = "bold")) 
dev.off()

dotplot(y, showCategory=50, #split=".sign",
        x = "NES",
        label_format = function(x) stringr::str_wrap(x, width=100),
        title = "Mog9 vs. OVA") + 
  #facet_grid(.~.sign) +
  #theme_prism(base_family = "Arial", base_size = 16) +
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  ) 

viewPathway(readable = T, foldChange = geneList)

gse <- gseGO(geneList=Treg_ranks, 
             ont ="BP", 
             keyType = "ENTREZID", 
             #nPerm = 10000, 
             minGSSize = 15, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.2, 
             verbose = TRUE, 
             by = "fgsea",
             OrgDb = org.Mm.eg.db, 
             #organism = "mouse",
             pAdjustMethod = "BH")

gse <- gseDGN(geneList=Treg_ranks, 
             pvalueCutoff = 0.2, 
             keyType = "ENTREZID",
             #verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             #organism = "mouse",
             pAdjustMethod = "BH")

check_gene_id(geneList, c(2020,4493,1030,3043,51115,54858))


dotplot(gse, showCategory=50, split=".sign",
        x = "NES",
        label_format = function(x) stringr::str_wrap(x, width=100),
        title = "Mog9 vs. OVA") + 
  facet_grid(.~.sign) +
  #theme_prism(base_family = "Arial", base_size = 16) +
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  ) 


#gse1 = gse
gse.df = as.data.frame(gse)
filter(gse, geneList == "Cyp2e1" )
gse@result$core_enrichment = list

class(gse1@result$core_enrichment)
 
heatplot(x, foldChange=geneList5, showCategory=25,# symbol = "dot",
         label_format = function(x) stringr::str_wrap(x, width=50)) +
  #theme_prism(base_family = "Arial", base_size = 16) +
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  ) 


plots <- VlnPlot(TemCD4, features = c("Ccl5", "Il7r", "Cd74", "Cd28", "Nkg7"
), 
#split.by = "sample",
group.by = "sample", pt.size = 0, combine = FALSE)

wrap_plots(plots, ncol = 2)
BiocManager::install("rWikiPathways")
library(rWikiPathways)
library(gson)
## downloaded from https://wikipathways-data.wmcloud.org/current/gmt/
gmt <- 'wikipathways-20230910-gmt-Homo_sapiens.gmt'
wp <- read.gmt.wp(downloadPathwayArchive(organism="Homo sapiens", format="gmt"))
data(DE_GSE8057)
colnames(DE_GSE8057)
colnames(wp)

length(TemCD4_D7_OvM[,1])
TemD7 = TemCD4_D7_OvM
TemD7$Day = "Day 7"
TemD9 = TemCD4_D9_OvM
TemD9$Day = "Day 9"
Tem = rbind(TemD7, TemD9)
Tem = Tem %>%
  distinct(gene,.keep_all = TRUE)
TemCD4ccgenes = data.frame("gene" = Tem$gene, "Day" = Tem$Day, "logFc" = Tem$avg_log2FC)


# Query AnnotationHub
ah <- AnnotationHub()
ah

unique(ah$species)
mouse_ens <- query(ah, c("Mus musculus", "EnsDb"))
mouse_ens <- mouse_ens[["AH64944"]]
annotations <- genes(mouse_ens, return.type = "data.frame")

colnames(EAEsig)[1] = "gene"
TemCD4ccgenes <- inner_join(TemCD4ccgenes, annotations, by=c("gene"="gene_name")) 
#TemCD4ccgenes = sort(TemCD4ccgenes, decreasing = TRUE)
#attach(TemCD4ccgenes)
#TemCD4ccgenes = TemCD4ccgenes[order(logFc, na.last = T),]
colnames(TemCD4ccgenes)
wp <- read.gmt.wp(downloadPathwayArchive(organism="Mus musculus", format="gmt"))
xx <- compareCluster(entrezid~Day, data=TemCD4ccgenes, fun = enricher,
                     TERM2GENE=wp[,c("wpid", "gene")], TERM2NAME=wp[,c("wpid", "name")], 
                     pvalueCutoff = 0.1)

TemCD4ccgenes$sign = "Up"
TemCD4ccgenes[which(TemCD4ccgenes$logFc < 0),]$sign = "Down"
xx <- compareCluster(entrezid~Day+sign, data=TemCD4ccgenes, fun = enricher,
                     TERM2GENE=wp[,c("wpid", "gene")], TERM2NAME=wp[,c("wpid", "name")], 
                     pvalueCutoff = 0.1)

png(file = "TemCD4_WPenrich_BMES.png",units = "in",width = 6, height = 4, res = 400)
clusterProfiler::dotplot(xx, x="Day", showCategory = 25,
              title = "CD4 Tem Cells\nWikiPathways Enrichment Analysis") +  facet_grid(~sign) +
  #aes(x=fct_relevel(time, c('0h', '2h', '6h', '24h'))) + 
  xlab(NULL) +
  scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y =element_text(size = 12, face = "bold")) 

dev.off()


TemD7 = TemCD8_D7_OvM
TemD7$Day = "Day 7"
TemD9 = TemCD8_D9_OvM
TemD9$Day = "Day 9"
Tem = rbind(TemD7, TemD9)
Tem = Tem %>%
  distinct(gene,.keep_all = TRUE)
TemCD8ccgenes = data.frame("gene" = Tem$gene, "Day" = Tem$Day, "logFc" = Tem$avg_log2FC)


# Query AnnotationHub
ah <- AnnotationHub()
ah

unique(ah$species)
mouse_ens <- query(ah, c("Mus musculus", "EnsDb"))
mouse_ens <- mouse_ens[["AH64944"]]
annotations <- genes(mouse_ens, return.type = "data.frame")

colnames(EAEsig)[1] = "gene"
TemCD8ccgenes <- inner_join(TemCD8ccgenes, annotations, by=c("gene"="gene_name")) 
#TemCD4ccgenes = sort(TemCD4ccgenes, decreasing = TRUE)
#attach(TemCD4ccgenes)
#TemCD4ccgenes = TemCD4ccgenes[order(logFc, na.last = T),]
colnames(TemCD8ccgenes)
wp <- read.gmt.wp(downloadPathwayArchive(organism="Mus musculus", format="gmt"))

TemCD8ccgenes$sign = "Up"
TemCD8ccgenes[which(TemCD8ccgenes$logFc < 0),]$sign = "Down"
xx <- clusterProfiler::compareCluster(entrezid~Day+sign, data=TemCD8ccgenes, fun = enricher,
                     TERM2GENE=wp[,c("wpid", "gene")], TERM2NAME=wp[,c("wpid", "name")], 
                     pvalueCutoff = 0.1)

png(file = "TemCD8_WPenrich_BMES.png",units = "in",width = 6, height = 4, res = 400)
clusterProfiler::dotplot(xx, x="Day",
        title = "Cytotoxic CD8 Tem Cells\nWikiPathways Enrichment Analysis") + facet_grid(~sign) +
  #aes(x=fct_relevel(time, c('0h', '2h', '6h', '24h'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y =element_text(size = 12, face = "bold")) 

dev.off()

TregD7 = Treg_D7_OvM
TregD7$Day = "Day 7"
TregD9 = Treg_D9_OvM
TregD9$Day = "Day 9"
Treg = rbind(TregD7, TregD9)
#remove duplicate rows
Treg = Treg %>%
  distinct(gene,.keep_all = TRUE)

Tregccgenes = data.frame("gene" = Treg$gene, "Day" = Treg$Day, "logFc" = Treg$avg_log2FC)


Tregccgenes <- inner_join(Tregccgenes, annotations, by=c("gene"="gene_name")) 
#TemCD4ccgenes = sort(TemCD4ccgenes, decreasing = TRUE)
#attach(TemCD4ccgenes)
#TemCD4ccgenes = TemCD4ccgenes[order(logFc, na.last = T),]
colnames(Tregccgenes)
#wp <- read.gmt.wp(downloadPathwayArchive(organism="Mus musculus", format="gmt"))

Tregccgenes$sign = "Up"
Tregccgenes[which(Tregccgenes$logFc < 0),]$sign = "Down"

head(Treg.genes.res)
Treg.genes.res1 =  

xx <- clusterProfiler::compareCluster(entrezid~Day+sign, data=Tregccgenes, fun = enricher,
                                      TERM2GENE=wp[,c("wpid", "gene")], TERM2NAME=wp[,c("wpid", "name")], 
                                      pvalueCutoff = 0.1)
getwd()
png(file = "Treg_WPenrich_updowninOVA.png",units = "in",width = 6, height = 4, res = 400)
clusterProfiler::dotplot(xx, x="Day",showCategory = 29,
                         title = "Treg Cells\nWikiPathways Enrichment Analysis") + 
  facet_grid(~sign) +
  #aes(x=fct_relevel(time, c('0h', '2h', '6h', '24h'))) + xlab(NULL) +
  scale_color_gradientn(colours=c("#b3eebe", "#46bac2", "#371ea3"),
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  guides(size = guide_legend(override.aes=list(shape=1))) +
  theme(panel.grid.major.y = element_line(linetype='dotted', color='#808080'),
        panel.grid.major.x = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y =element_text(size = 12, face = "bold")) 

dev.off()

pp

gse <- gseGO(geneList=geneList1, 
             ont ="BP", 
             keyType = "SYMBOL", 
             #nPerm = 10000, 
             minGSSize = 10, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")
dotplot(gse, showCategory=50, split=".sign",
        x = "GeneRatio",
        label_format = function(x) stringr::str_wrap(x, width=100),
        title = "Mog9 vs. OVA") + 
  facet_grid(.~.sign) +
  #theme_prism(base_family = "Arial", base_size = 16) +
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  ) 

inputList <- list(Day7 = geneList1, Day9 = geneList2)

inputList <- list(Day7 = geneList3, Day9 = geneList4)
test.out <- compareCluster(geneClusters=inputList,  fun = "gseGO", 
                           ont ="BP", 
                           keyType = "ENTREZID", 
                           #nPerm = 10000, 
                           minGSSize = 10, 
                           maxGSSize = 1000, 
                           pvalueCutoff = 0.1, 
                           verbose = TRUE, 
                           OrgDb = org.Mm.eg.db, 
                           pAdjustMethod = "BH")


test.out <- compareCluster(geneClusters=inputList,  fun = "gseWP", 
                           organism = "Mus musculus", pvalueCutoff = 0.1)

gseW=gseWP(geneList2, organism = "Mus musculus", pvalueCutoff = 0.1)
gseK = gseKEGG(geneList3, organism = "mmu", pvalueCutoff = 0.1)

dotplot(gseK, showCategory=100, split=".sign") + facet_grid(.~.sign)
clusterProfiler::dotplot(test.out, showCategory=50, split=".sign",
        #x = "GeneRatio",
        label_format = function(x) stringr::str_wrap(x, width=100),
        title = "Mog9 vs. OVA") + 
  facet_grid(.~.sign) +
  #theme_prism(base_family = "Arial", base_size = 16) +
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  ) 
gse.df = as.data.frame(test.out)

#Rename Idents in larger seurat object ----
# Generate a new column called sub_cluster in the metadata
TNBC$sub_cluster <- as.character(Idents(TNBC))

# Change the information of cells containing sub-cluster information
Idents(TsubClusters)
TNBC$sub_cluster[Cells(TsubClusters)] <- as.character(Idents(TsubClusters))
DimPlot(TNBC, group.by = "sub_cluster", label = T, repel = T)
saveRDS(TNBC, "LMR_16renamedclusters_wTcellsubcluster.rds")



# Old Code ----
# T cell data ----
Idents(TNBC)
Tcell_data = subset(TNBC,idents = c("CD4 T", "CD8 T", "TH17"))
DGE_cell_selection <- FindAllMarkers(Tcell_data, log2FC.threshold = 0.2, test.use = "wilcox",
                                     min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                     assay = "RNA")

DGE_cell_selection %>%
  group_by(cluster) %>%
  top_n(-5, p_val) -> top5_cell_selection

DGE_cell_selection %>%
  group_by(cluster) %>%
  top_n(-25, p_val) -> top25_cell_selection

top25_cell_selection = top25_cell_selection[top25_cell_selection$p_val_adj <= 0.05,]

RidgePlot(Tcell_data, features = as.character(unique(top25_cell_selection$gene))[1:6], ncol = 2)

png(file = "topDETcell_poster.png",
    units = "in",width = 45, height = 12, res = 400)
VlnPlot(Tcell_data, features = as.character(unique(top25_cell_selection$gene))[1:6], ncol = 3,
          split.by = "sample")  & theme_prism(base_size = 34, base_family = "Arial") & NoLegend() 
dev.off()

png(file = "topDETcellLegend_poster.png",
    units = "in",width = 45, height = 20, res = 400)
VlnPlot(Tcell_data, features = as.character(unique(top25_cell_selection$gene))[1], ncol = 2,
        split.by = "sample")  & theme(legend.position = "top")
dev.off()

VlnPlot(Tcell_data, features = as.character(unique(top25_cell_selection$gene)),
        ncol = 5, group.by = "sample", assay = "SCT", pt.size = 0.1)

png(file = "tcellclusters_sctype.png",
    units = "in",width = 40, height = 20, res = 400)
VlnPlot(Tcell_data, features = c("Cd3e", "Cd4", "Cd8a", "Itga4", "Vcam1","Itgal", "Itgb1", "Ccr7"),
        ncol = 5, assay = "RNA", pt.size = 0.1, log = T, split.by="sample" )
dev.off()

png(file = "tcellclusters_sctype.png",
    units = "in",width = 40, height = 20, res = 400)
RidgePlot(Tcell_data, features = c("Cd3e", "Cd4", "Cd8a", "Il2ra", "Vcam1","Itgal", "Itgb1", "Ccr7"),
       assay = "RNA")
dev.off()

png(file = "tcellclusters_poster.png",
    units = "in",width = 60, height = 10, res = 400)
VlnPlot(Tcell_data, features = c("Cd3e", "Cd4", "Cd8a", "Ccr7", "Ifngr1","Ctla4"),
          assay = "RNA", split.by = "Antigen", ncol = 6) & theme_prism(base_size = 34, base_family = "Arial") & NoLegend() 
dev.off()

png(file = "tcellclustersLegend_poster.png",
    units = "in",width = 30, height = 20, res = 400)
VlnPlot(Tcell_data, features = c("Cd3e"),
        assay = "RNA", split.by = "Antigen") + theme_prism(base_size = 34, base_family = "Arial") +
  theme(legend.position = "top")
dev.off()

Idents(Tcell_data) = Tcell_data[["sample"]]

MOG7v9 <- FindMarkers(Tcell_data, ident.1 = "D7 MOG", ident.2 = "D9 MOG")
MOG7v9p <- FindMarkers(Tcell_data, ident.1 = "D7 MOG", ident.2 = "D9 MOG", only.pos = T)
setdiff(rownames(MOG7v9), rownames(MOG7v9p))
setdiff(rownames(MOG7v9p), rownames(MOG7v9))

CD4Tcells = subset(TNBC,idents = c("CD4 T"))
CD8Tcells = subset(TNBC,idents = c("CD8 T"))

Idents(CD4Tcells) = CD4Tcells[["sample"]]
Idents(CD8Tcells) = CD8Tcells[["sample"]]
CD4MOG7v9 <- FindMarkers(CD4Tcells, ident.1 = "D7 MOG", ident.2 = "D9 MOG")
#MOG7v9p <- FindMarkers(Tcell_data, ident.1 = "D7 MOG", ident.2 = "D9 MOG", only.pos = T)

CD4OVA7v9 <- FindMarkers(CD4Tcells, ident.1 = "D7 OVA", ident.2 = "D9 OVA")

RelevantCD4genes = setdiff(rownames(CD4MOG7v9), rownames(CD4OVA7v9))


VlnPlot(Tcell_data, features = "Ccr7", split.by = "sample")

CD4MOG7v9 = CD4MOG7v9[RelevantCD4genes,]
CD4MOG7v9 = CD4MOG7v9[CD4MOG7v9$avg_log2FC > 0.25,]
CD4MOG7v9 = CD4MOG7v9[complete.cases(CD4MOG7v9), ]
CD4MOG7v9 = CD4MOG7v9[!grepl("^Rp|^Rps|^H", rownames(CD4MOG7v9)),]
CD4MOG7v9 = CD4MOG7v9[which(CD4MOG7v9$avg_log2FC>0.65),]

png(file = "tcellclusters_poster.png",
    units = "in",width = 60, height = 10, res = 400)
VlnPlot(CD4Tcells, features = rownames(CD4MOG7v9),
        assay = "RNA", split.by = "sample") & theme_prism(base_size = 14, base_family = "Arial") & NoLegend() 
dev.off()

png(file = "tcellclustersLegend_poster.png",
    units = "in",width = 30, height = 20, res = 400)
VlnPlot(Tcell_data, features = c("Cd3e"),
        assay = "RNA", split.by = "Antigen") + theme_prism(base_size = 34, base_family = "Arial") +
  theme(legend.position = "top")
dev.off()

#CD8 ----
CD8MOG7v9 <- FindMarkers(CD8Tcells, ident.1 = "D7 MOG", ident.2 = "D9 MOG")
#MOG7v9p <- FindMarkers(Tcell_data, ident.1 = "D7 MOG", ident.2 = "D9 MOG", only.pos = T)

CD8OVA7v9 <- FindMarkers(CD8Tcells, ident.1 = "D7 OVA", ident.2 = "D9 OVA")

RelevantCD8genes = setdiff(rownames(CD8MOG7v9), rownames(CD8OVA7v9))


CD8MOG7v9 = CD8MOG7v9[RelevantCD8genes,]
CD8MOG7v9 = CD8MOG7v9[abs(CD8MOG7v9$avg_log2FC)>0.25,]
CD8MOG7v9 = CD8MOG7v9[complete.cases(CD8MOG7v9), ]
CD8MOG7v9 = CD8MOG7v9[!grepl("^Rp|^Rps|^H", rownames(CD8MOG7v9)),]
attach(CD8MOG7v9)
CD8MOG7v9$avg_log2FC = abs(CD8MOG7v9$avg_log2FC)
CD8MOG7v9 = CD8MOG7v9[abs(CD8MOG7v9$avg_log2FC)>0.65,]
CD8MOG7v9 = CD8MOG7v9[CD8MOG7v9$p_val_adj<0.05,]

png(file = "topDECD8Tcell_poster.png",
    units = "in",width = 60, height = 10, res = 400)
VlnPlot(CD8Tcells, features = rownames(CD8MOG7v9), ncol = 6,
        split.by = "sample")  & theme_prism(base_size = 34, base_family = "Arial") & NoLegend() 
dev.off()

png(file = "topDETcellLegend_poster.png",
    units = "in",width = 45, height = 20, res = 400)
VlnPlot(CD8Tcells, features = rownames(CD8MOG7v9)[1], ncol = 2,
        split.by = "sample")  & theme(legend.position = "top")
dev.off()

#Itga4: VLA4 subunit  CD49d (alpha 4) 
#Itgb1: VLA4 subunit CD29 

#Vcam1: On endothelial cells and binds to VLA4

#Itga1: subunit Lfa1 leukocyte adhesion molecule to cns

png(file = "CD4Tcells_VLA4.png",
    units = "in",width = 25, height = 30, res = 400)
VlnPlot(CD4Tcells, features = c("Cd4","Cd8a" ,"Itga4","Itgb1","Itgal", "Itgb2" ),
        ncol = 2, group.by = "sample", pt.size = 0.1) & 
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=30,face="bold"),
        axis.text.y = element_text(size=30,face="bold"),
        axis.title.y = element_text(size=30,face="bold"))
dev.off()

png(file = "CD8Tcells_VLA4.png",
    units = "in",width = 25, height = 30, res = 400)
VlnPlot(CD8Tcells, features = c("Cd4","Cd8a" ,"Itga4","Itgb1","Itgal", "Itgb2" ),
        ncol = 2, group.by = "sample", assay = "SCT", pt.size = 0.1)  & 
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=30,face="bold"),
        axis.text.y = element_text(size=30,face="bold"),
        axis.title.y = element_text(size=30,face="bold"))
dev.off()

VlnPlot(data.filt, features = "Vcam1")

CD4D7MvO <- FindMarkers(CD4Tcells, ident.1 = "D7 MOG", ident.2 = "D7 OVA")
CD4D9MvO <- FindMarkers(CD4Tcells, ident.1 = "D9 MOG", ident.2 = "D9 OVA")

CD4D7MvO = CD4D7MvO[CD4D7MvO$p_val_adj <= 0.05,]
CD4D9MvO = CD4D9MvO[CD4D9MvO$p_val_adj <= 0.05,]

png(file = "CD4_D7_MvO.png",
    units = "in",width = 20, height = 5, res = 400)
VlnPlot(CD4Tcells, features = rownames(CD4D7MvO),
        ncol = 6, group.by = "sample", assay = "RNA", pt.size = 0.1)
dev.off()

VlnPlot(data.filt, features = rownames(CD4D7MvO),
        ncol = 3, assay = "SCT", pt.size = 0.1)

png(file = "CD4_D9_MvO.png",
    units = "in",width = 20, height = 10, res = 400)
VlnPlot(CD4Tcells, features = rownames(CD4D9MvO),
        ncol = 5, group.by = "sample", assay = "SCT", pt.size = 0.1)
dev.off()


#Gene Set Analysis ----
# Load additional packages
library(enrichR)

# Check available databases to perform enrichment (then choose one)
enrichR::listEnrichrDbs()
new.cluster.ids

#subset of innate cells only
innateCells = subset(TNBC,idents = c("Monocytes","Macrophage","Complement Macrophage","Neutrophil","Dendritic Cell",
                                     "Chemoattractant Monocyte","Immature DC"))

innateCells <- AddMetaData(
  object = innateCells,
  metadata = Idents(innateCells),
  col.name = 'cellType'
)



# differentially expressed genes from innate subset between Day 7 Mog and Day 7 OVA
Idents(innateCells) = innateCells[["sample"]]
innate7MvO<- FindMarkers(innateCells, ident.1 = "D7 MOG", ident.2 = "D7 OVA") 


innate7MvO = innate7MvO[abs(innate7MvO$p_val_adj)<0.05,]
innate7MvO = innate7MvO[complete.cases(innate7MvO), ]
#remove ribosomal and histone genes
innate7MvO = innate7MvO[!grepl("^Rp|^Rps|^H", rownames(innate7MvO)),]


Idents(innateCells) = innateCells[["cellType"]]

VlnPlot(innateCells, features = rownames(innate7MvO)[1:6],
        assay = "RNA") & theme_prism(base_size = 14, base_family = "Arial") & NoLegend()

# Perform enrichment

#c("BioPlanet_2019","Reactome_2022") Immune relevant gene sets
enrich_resultsInnate <- enrichr(genes =  c(rownames(innate7MvO)), databases = "BioPlanet_2019")
enrich_resultsInnate = as.data.frame(enrich_resultsInnate[[1]])
#enrich_results = enrich_results[-which(enrich_results$Term == "Gastrin pathway"),] #remove Gastrin Pathway

png(file = "InnateCellsDay7enrichment.png",
    units = "in",width = 32, height = 13, res = 400)
plotEnrich(enrich_resultsInnate, showTerms = 20, numChar = 80, y = "Count", orderBy = "Adjusted.P.value") & 
  ggtitle("Functional Enrichment of Innate Cells at Day 7") &
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=20,face="bold"),
        axis.text.y = element_text(size=40,face="bold"),
        axis.title.y = element_text(size=20,face="bold"),
        legend.text = element_text(size=30,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=30,face="bold",
                    
                                   ))
dev.off()

DCs = subset(TNBC,idents = c("Dendritic Cell","Immature DC"))

DCs <- AddMetaData(
  object = DCs,
  metadata = Idents(DCs),
  col.name = 'cellType'
)



# differentially expressed genes from innate subset between Day 7 Mog and Day 7 OVA
Idents(DCs) = DCs[["sample"]]
DC7MvO<- FindMarkers(DCs, ident.1 = "D7 MOG", ident.2 = "D7 OVA") 


#DC7MvO = DC7MvO[abs(DC7MvO$p_val_adj)<0.05,]
DC7MvO = DC7MvO[complete.cases(DC7MvO), ]
#remove ribosomal and histone genes
DC7MvO = DC7MvO[!grepl("^Rp|^Rps|^H", rownames(DC7MvO)),]

#top DEG based on LFC
DC7MvOtop = DC7MvO
DC7MvOtop$avg_log2FC = abs(DC7MvOtop$avg_log2FC)
DC7MvOtop = DC7MvOtop[abs(DC7MvOtop$avg_log2FC)>0.65,]


Idents(DCs) = DCs[["cellType"]]

VlnPlot(DCs, features = rownames(DC7MvOtop)[1:6],
        assay = "RNA") & theme_prism(base_size = 14, base_family = "Arial") & NoLegend()
VlnPlot(DCs, features = rownames(DC7MvOtop)[1:6],
        assay = "RNA", split.by = "sample") & theme_prism(base_size = 14, base_family = "Arial") 

# Perform enrichment

#c("BioPlanet_2019","Reactome_2022") Immune relevant gene sets
enrich_resultsDCs <- enrichr(genes =  c(rownames(DC7MvO)), databases = "BioPlanet_2019")
enrich_resultsDCs = as.data.frame(enrich_resultsDCs[[1]])
#enrich_results = enrich_results[-which(enrich_results$Term == "Gastrin pathway"),] #remove Gastrin Pathway

png(file = "DCsDay7enrichment.png",
    units = "in",width = 32, height = 13, res = 400)
plotEnrich(enrich_resultsDCs, showTerms = 20, numChar = 80, y = "Count", orderBy = "Adjusted.P.value") & 
  ggtitle("Functional Enrichment of DCs at Day 7") &
  theme(text = element_text(face = "bold"),
        plot.title = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=20,face="bold"),
        axis.text.y = element_text(size=40,face="bold"),
        axis.title.y = element_text(size=20,face="bold"),
        legend.text = element_text(size=30,face="bold"),
        legend.key.size = unit(1, 'cm'),
        legend.title = element_text(size=30,face="bold",
                                    
        ))
dev.off()


#T cell Clustering IDs ------
Idents(Tcell_data) = Tcell_data[["sample"]]
Tcell_data<-FindVariableFeatures(Tcell_data)
s.genes <- str_to_sentence(cc.genes$s.genes)
g2m.genes <- str_to_sentence(cc.genes$g2m.genes)
Tcell_data <- CellCycleScoring(Tcell_data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Tcell_data <- ScaleData(Tcell_data, vars.to.regress=c("S.Score","G2M.Score","nUMI","percent.mt"), verbose = TRUE)
Tcell_data<-FindVariableFeatures(Tcell_data)
Tcell_data <- RunPCA(Tcell_data, npcs = 30, verbose = TRUE)
# t-SNE and Clustering
Tcell_data <- RunUMAP(Tcell_data, reduction = "pca", dims = 1:20)
Tcell_data <- FindNeighbors(Tcell_data, reduction = "pca", dims = 1:20)
Tcell_data<-FindClusters(Tcell_data,resolution=.3)


DimPlot(Tcell_data, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T)

png(file = "UMAPTcells.png",
    units = "in",width = 16, height = 11, res = 400)
DimPlot(Tcell_data, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T) + 
  theme_prism(base_family = "Arial", base_size = 24) + theme(legend.text = element_text(size = 28))
#+ NoLegend()
dev.off()


Tcell.markers<-FindAllMarkers(Tcell_data,verbose = TRUE,only.pos = TRUE)



dim(Tcell.markers)
table(Tcell.markers$cluster)
table(Tcell_data$seurat_clusters)

# Frequency plots
Idents(Tcell_data) = Tcell_data[["seurat_clusters"]]
TCellTypeTable = as.data.frame(proportions(table(Idents(Tcell_data),Tcell_data@meta.data$sample), margin = 2))
TCellTypeTable$Freq = TCellTypeTable$Freq*100
levels(TCellTypeTable$Var1)


getPalette = colorRampPalette(brewer.pal(9, "Set1"))


png(file = "TCellProportions.png",
    units = "in",width = 9, height = 5, res = 400)
ggplot(TCellTypeTable,aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism(base_family = "Arial", base_size = 16)+
  labs(x = "Condition", y= "Percent of Cells",title = "T Cell Subset Frequencies Across Samples") + 
  scale_fill_manual(values = getPalette(8), name = "Cell Type") + 
  theme(legend.text = element_text(size = 16))
dev.off()


top10_markers <- as.data.frame(Tcell.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC))
top10_markers

VlnPlot(Tcell_data, features = as.character(unique(top3_markers$gene)), ncol = 10, group.by = "seurat_clusters",
        assay = "RNA", pt.size = 0)



DotPlot(Tcell_data, features = c("Cd3e", "Cd4", "Cd8a", 
                                  "Cd19", "Ighm", 
                                 "Ly6c2", "Ccr2", "Ccr7",
                                 "Ccr5", "Ctla4", "Lag3",
                                "Izumo1r", #Fr-4
                                "Nt5e", #CD73
                                "Foxp3", "Gata3", "Il17a",
                                "Itgb1", "Ccl5",
                                "Il2ra", "Ifng","Ifngr1", "Pdcd1",
                                "Rorc", "Tbx21", #Th17
                                "Havcr2", "Sell", #CD62L
                                "Klrg1", #CD8 memory
                                "Il10", "Il10ra",
                                "Slamf6", "Xcl1", "Gzmb", "Gzmk",
                                "Cx3cr1", "Cxcr5","Tox"
                                
)) + 
  RotatedAxis()


Tcellmarkers = c(
  "Cd3e", "Cd8a", "Cd4",
  #naive T cells:
  "Ccr7", "Cd28",
  #cytotoxic CD8:
  "Prf1","Gzmb",
  #Th1:
  "Ccr4","Cxcr3", "Ccr5", "Il12rb1","Il12rb2", "Ifng","Ifngr1", "Tbx21", #Tbet 
  #Th2:
  "Il4","Il4ra", "Il5","Il5ra", "Gata3",
  #Th9:
  "Ccr3", "Ccr6", "Spi1","Il9r",#"Il22ra1",
  #Th17: Ccr6, Ccr4 Nk1.1
  "Klrb1c", "Il17a","Il17f","Il6ra","Tgfb1","Il23a",
  #Treg: CD127 = Il7r
  "Il7r", "Il2ra", "Ctla4", "Il2", "Foxp3", "Il10","Il10ra",
  #Tfh:
  "Cxcr5", "Cd40lg", "Icos", "Bcl6", "Il21r",
  #gammadelta T cell:
  "Il23r", "Rorc","Rora",
  #Memory:
  "Cd44","Sell", #CD62L
  "Cd27", "Itgae", #Cd103
  "Lag3",
  "Izumo1r", #Fr-4
  "Nt5e", #CD73
  "Slamf6", "Xcl1","Cx3cr1", "Tox", "Gzmk", "Pdcd1",
  "Trgv2", "Trdc", "Trdv4"
)
png(file = "Tcell_clusterbyhanddotplot.png",
    units = "in",width = 20, height = 10, res = 400)
DotPlot(Tcell_data, features = Tcellmarkers#, cols = c("blue", "red")
        ) + 
  RotatedAxis()
dev.off()
VlnPlot(Tcell_data, features = c("Foxp3","Itgb1", "Cd8a"), split.by = "sample")


# Cluster 0 is Naive CD4
# Cluster 1 is CD4+ Effector Memory Cd44+ CD62L- (Sell-)
# Cluster 2 is Naive CD8
# Cluster 3 is Cytotoxic CD8
# Cluster 4 is Treg
# Cluster 5 is Double Positive
# Cluster 6 are gamma/delta 
# Cluster 7 are Double Positive (2)

Tnew.cluster.ids<-c("Naive CD4","CD4+ Effector Memory","Naive CD8","Cytotoxic CD8","Treg","Double Positive","gd"," Double Positive (2)")
names(Tnew.cluster.ids) <- levels(Tcell_data)
Tcell_data <- RenameIdents(Tcell_data, Tnew.cluster.ids)

DGE_Tcell_selection <- FindAllMarkers(Tcell_data, log2FC.threshold = 0.2, test.use = "wilcox",
                                     min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50,
                                     assay = "RNA")

DGE_Tcell_selection %>%
  group_by(cluster) %>%
  top_n(-5, p_val) -> top5_Tcell_selection

Idents(Tcell_data) = Tcell_data[["sample"]]

MOG7v9_T <- FindMarkers(Tcell_data, ident.1 = "D7 MOG", ident.2 = "D9 MOG")

MOG7v9_T = MOG7v9_T[MOG7v9_T$avg_log2FC > 0.25,]
MOG7v9_T = MOG7v9_T[complete.cases(MOG7v9_T), ]
MOG7v9_T = MOG7v9_T[!grepl("^Rp|^Rps|^H", rownames(CD4MOG7v9)),]
MOG7v9_T = MOG7v9_T[which(MOG7v9_T$avg_log2FC>0.65),]

VlnPlot(Tcell_data, features = rownames(MOG7v9_T)[1:6], split.by = "sample")

head(Tcell_data[[]])
Tcell_data[["CellTypes"]] = Idents(Tcell_data)
Idents(Tcell_data) = Tcell_data[["CellTypes"]]


TregDE = FindMarkers(Tcell_data, ident.1 = "D7 MOG", ident.2 = "D9 MOG",
                     group.by = "sample", subset.ident = "Cytotoxic CD8")


#Functional Enrichment
library(AnnotationHub)
library(ensembldb)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(genekitr)

ah <- AnnotationHub()
ah
library(msigdbr)
library(enrichR)

# Query AnnotationHub
mouse_ens <- query(ah, c("Mus musculus", "EnsDb"))
mouse_ens <- mouse_ens[["AH64944"]]
annotations <- genes(mouse_ens, return.type = "data.frame")


Sig = as.data.frame(rownames(TregDE[which(TregDE$p_val_adj < 0.1),]))

colnames(Sig) = "gene"

res_ids <- inner_join(Sig, annotations, by=c("gene"="gene_name")) 

msigdbr_species()
m_df <- msigdbr(species = "Mus musculus")
head(m_df, 2) %>% as.data.frame
m_t2g <- msigdbr(species = "Mus musculus", category = "C2") %>% 
  dplyr::select(gs_name, entrez_gene)
head(m_t2g)

em <- enricher(res_ids$entrezid, TERM2GENE=m_t2g)
head(em)
dotplot(em, showCategory = 29, title = "GO Enrichment Analysis of Spleen DEGs between PBS and NP")
barplot(em, showCategory=29, 
        title = "GO Enrichment Analysis of Spleen DEGs between PBS and NP",
        label_format = function(x) stringr::str_wrap(x, width=60),
        font.size = 10) 

geneList = TregDE$avg_log2FC
names(geneList) = as.character(unique(res_ids$entrezid))

## feature 3: decreasing orde
geneList = sort(geneList, decreasing = TRUE)

em2 <- clusterProfiler::GSEA(geneList, TERM2GENE = m_t2g)
head(em2)

gsea_import = importCP(em2,type = 'other')

plotGSEA(em2, plot_type = "classic", show_pathway = pathways, show_gene = genes)

dotplot(em2, showCategory=29, 
        title = "GO Enrichment Analysis of Spleen DEGs between PBS and NP",
        label_format = function(x) stringr::str_wrap(x, width=60),
        font.size = 10) 


ego <- enrichGO(gene = as.character(res_ids$gene_id), 
                ont = "BP",
                #universe = as.character(allgenes$gene_id),
                keyType = "ENSEMBL",
                OrgDb = org.Mm.eg.db, 
                pvalueCutoff = 0.05, 
                readable = TRUE)

tab.go <- data.frame(ego)


png(file = paste(path, "GO_PP_DEGsPBSvNP_Genes_dotplot.png",sep=""),
    units = "in",width = 6, height = 12, res = 400)
dotplot(ego, showCategory = 29, title = "GO Enrichment Analysis of PP DEGs between PBS and NP")
dev.off()

png(file = paste(path, "GO_barplot_PP_DEGsPBSvNP.png",sep=""),
    units = "in",width = 10, height = 4, res = 400)
barplot(ego, showCategory=29, 
        title = "GO Enrichment Analysis of PP DEGs between PBS and NP",
        label_format = function(x) stringr::str_wrap(x, width=60),
        font.size = 10) 
#ggsave(file = paste(path, "GO_barplot_Mv09t_uni.svg",sep=""), width = 20, height = 10, pointsize = 12)
dev.off()

## Add similarity matrix to the termsim slot of enrichment result
ego <- enrichplot::pairwise_termsim(ego)

png(file = paste(path, "GO_emapplot_M14_9t.png",sep=""),
    units = "in",width = 12, height = 12, res = 400)
## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 30, font = 3, repel = T,
         title = "GO Term Clustering MOG D14 vs D9")
dev.off()

deGenes <- unlist(mget(as.character(res_ids$gene_id), envir=org.Mm.egENSEMBL2EG,
                       ifnotfound = NA))

#write.csv(res_ids$gene,paste(path,"/EAE_BI_cohort1DEmogTvC.csv", sep = ""), row.names = T)

#geneUniverse <- unlist(mget(as.character(allgenes$gene_id), envir=org.Mm.egENSEMBL2EG,ifnotfound = NA))

e.kegg = enrichKEGG(gene = deGenes, 
                    organism = 'mmu',
                    #universe = geneUniverse,
                    pvalueCutoff = 0.05)
tab.kegg <- as.data.frame(e.kegg)

png(file = paste(path, "KEGG_dotplot.png",sep=""),
    units = "in",width = 6, height = 12, res = 400)
dotplot(e.kegg, showCategory = 12, title = "KEGG Enrichment Analysis")
dev.off()

# Perform enrichment

#c("BioPlanet_2019","Reactome_2022", "Elsevier_Pathway_Collection","GO_Biological_Process_2023") Immune relevant gene sets
enrich_results <- enrichr(genes =  Sig$gene, databases = "BioPlanet_2019")

enrich_results= as.data.frame(enrich_results[[1]])
#enrich_results = enrich_results[-which(enrich_results$Term == "Gastrin pathway"),] #remove Gastrin Pathway

png(file = "enrichr_barplot_PP_DEGsPBSvNP.png",
    units = "in",width = 32, height = 13, res = 400)
plotEnrich(enrich_results, showTerms = 29, numChar = 80, y = "Count", orderBy = "Adjusted.P.value") & 
  ggtitle("Functional Enrichment of Innate Cells at Day 7")

dev.off()

# automated cell type annotation azimuth ----
#devtools::install_github("satijalab/azimuth")

library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

pbmcsca <- RunAzimuth(TNBC, reference = "pbmcref")
pbmcsca[[]]
p1 <- DimPlot(pbmcsca, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(pbmcsca, group.by = "predicted.celltype.l3", label = TRUE, label.size = 3) + NoLegend()
p1 + p2

unique(pbmcsca$predicted.celltype.l2)

Idents(pbmcsca) = "predicted.celltype.l2"

p1 <- FeaturePlot(pbmcsca, features = "CCR7")
p2 <- FeaturePlot(pbmcsca, features = "FCGR3A")
p3 <- VlnPlot(pbmcsca, features = "CD8A", group.by = "predicted.celltype.l2",split.by = "sample")
p4 <- VlnPlot(pbmcsca, features = "CD4", group.by = "predicted.celltype.l2",split.by = "sample")

p1 + p2 + p3 + p4 + plot_layout(ncol=1)



#joe's other code ----

########NICHES is what we're going to try now#############


library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(NICHES)
library(dplyr)
library(ggplot2)


## Convert to Cell-Cell Signaling Atlas
# Here, we use the human-specific OmniPath database as our ground-truth ligand-receptor mechanism database, which includes some mechanisms with more than one subunit on both the sending and receiving arms:


TNBC<-subset(TNBC,idents=c("0","2","11","9","7"))

Idents(TNBC)<-c("Monocytes","CD8 T Cell","Dendritic Cell","CD4 T Cell","B Cell")

scc <- RunNICHES(TNBC,
                 assay = 'RNA',
                 species = 'mouse',
                 LR.database = 'omnipath',
                 CellToCell = T)


## Visualize Cell-Cell Signaling Relationships using UMAP
# We can visualize cell-cell signaling relationships quickly using tSNE or UMAP as follows:

demo <- scc$CellToCell
demo <- ScaleData(demo)
demo <- FindVariableFeatures(demo)
demo <- RunPCA(demo)
ElbowPlot(demo,ndims=50)
PCHeatmap(demo,dims = 1:6,balanced = T,cells = 100)
demo <- RunUMAP(demo,dims = 1:6)
DimPlot(demo,reduction = 'umap',group.by = 'VectorType',label = F)

# This allows us to see some patterns and segregation of cell-cell relationships, however, this doesn't look particularly clean. This is because the cell types that we have chosen here overlap significantly in the mechanisms that they employ. If we want to further clarify patterns that might not otherwise be obvious, we can run NICHES on the same data but imputed to amplify low-strength signal and remove gene-sampling-related clustering artifacts:
  
 
imputed <- SeuratWrappers::RunALRA(TNBC)
scc.imputed <- RunNICHES(imputed,
                         assay = 'alra',
                         species = 'mouse',
                         LR.database = 'omnipath',
                         CellToCell = T)
demo.2 <- scc.imputed$CellToCell
demo.2 <- ScaleData(demo.2)
demo.2 <- FindVariableFeatures(demo.2)
demo.2 <- RunPCA(demo.2)
ElbowPlot(demo.2,ndims=50)
PCHeatmap(demo.2,dims = 1:6,balanced = T,cells = 100)
demo.2 <- RunUMAP(demo.2,dims = 1:30)
demo.2<-NormalizeData(demo.2)
demo.2<-FindNeighbors(demo.2)
demo.2<-FindClusters(demo.2)
DimPlot(demo.2,reduction = 'umap',group.by = 'VectorType',label = F) +NoLegend()

#This looks cleaner, removes sampling-related artifacts, and allows clearer observation of relationships which occupy their own region of mechanism-space versus overlap with other cell-cell relationships. Differential expression can be run to identify prominent marker mechanisms as follows:
  
# Find markers
mark <- FindAllMarkers(demo.2)
GOI_niche <- mark %>% group_by(cluster) %>% top_n(5,p_val)
DoHeatmap(demo.2,features = unique(GOI_niche$gene)) 
  scale_fill_gradientn(colors = c("grey","white", "blue"))























ll<-TNBC[TNBC]




Cell_Populations<-table(Idents(TNBC),TNBC$sample)

plot(Cell_Populations,color=c("red","orange","yellow","green"))









pbmc<-readRDS("TNBC_scRNAseq_AZIZI.rds")

df2regulon <- function(df) {
  regulon_list = split(df, df$tf)
  
  viper_regulons = lapply(regulon_list, function(regulon) {
    tfmode = stats::setNames(regulon$mor, regulon$target)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })
  
  return(viper_regulons)
}







pbmc<-lung



library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)

dorothea_regulon_human <- get(data("dorothea_mm", package = "dorothea"))

regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C","D"))


Idents(pbmc) = as.factor(melanoma_cell_types)


pbmc<-lung#readRDS(file="/home/jtdecker/TNBC_scRNAseq.rds")

pbmc<-TNBC

DefaultAssay(pbmc)<-"RNA"

pbmc<-FindVariableFeatures(pbmc)
pbmc<-TNBC

pbmc<- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

pbmc <- ScaleData(pbmc, vars.to.regress=c("S.Score","G2M.Score"), verbose = TRUE)
pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)

UMAPPlot(pbmc)

pbmc <- run_viper(pbmc, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 1, 
                                 verbose = TRUE))


viper_scores_df <- GetAssayData(TNBC,  
                                assay = "dorothea") %>%
  data.frame() %>%
  t()


viper_scores_df<-viper_scores_df[which(Idents(pbmc)=="T Cells"),]

Timepoint<-as.factor(pbmc$Time[which(Idents(pbmc)=="T Cells")])

## We create a data frame containing the cells and their clusters
CellsClusters <- data.frame(cell = row.names(viper_scores_df),
                            cell_type = Timepoint,
                            stringsAsFactors = FALSE)


## We create a data frame with the Viper score per cell and its clusters
viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)

## We summarize the Viper scores by cellpopulation
summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  plotly::summarise(avg = mean(activity),
            std = sd(activity))




highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  plotly::mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(20*4, var) %>%
  distinct(tf)



## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=14, 
                       fontsize_row = 10, 
                       color=my_color, breaks = my_breaks, 
                       main = "Lung T Cell Variable TFs", 
                       treeheight_col = 0,  border_color = NA,cluster_cols = F) 




pbmc<-subset(immune.combined,subset=seurat_clusters%in%c("9","8","2"))


ifnb.list <- SplitObject(pbmc, split.by = "sample")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:10,k.filter=1,k.score=1,k.anchor=1)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20,k.weight=20)
DefaultAssay(immune.combined) <- "integrated"

immune.combined<-pbmc


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

immune.combined <- ScaleData(immune.combined, vars.to.regress=c("S.Score","G2M.Score","nUMI","percent.mt"), verbose = TRUE)
immune.combined<-FindVariableFeatures(immune.combined)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = TRUE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined<-FindClusters(immune.combined,resolution=.3)

DefaultAssay(immune.combined)="RNA"


pbmc<-immune.combined

DefaultAssay(object = pbmc) <- "integrated"
pbmc <- ScaleData(pbmc,verbose=TRUE)#, vars.to.regress=c("S.Score","G2M.Score","nUMI","percent.mt"), verbose = TRUE)

pbmc <- RunPCA(pbmc, features = rownames(pbmc), verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:20, verbose = FALSE)
pbmc <- FindClusters(pbmc,resolution=.2)# verbose = FALSE)

pbmc <- RunUMAP(pbmc, dims = 1:10, umap.method = "uwot", metric = "cosine")

DefaultAssay(object=pbmc)<-"RNA"

pbmc<-NormalizeData(pbmc)

Tcells[]

New.idents<-c("Naive CD4","Exhausted CD8","Treg","Fibroblast H","Nk","Cytotoxic CD8","Activated CD4","Endothelial H","B H","Mast H")


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, 
                               logfc.threshold = 0.25, verbose = TRUE)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


library(cowplot)
p1 <- DimPlot(TNBC, reduction = "umap", group.by = "sample")
p2 <- DimPlot(TNBC, reduction = "umap", label = TRUE)
plot_grid(p1, p2)


saveRDS(pbmc,file="/home/jtdecker/TNBC_scRNAseq.rds")

MO<-readRDS(file="/home/jtdecker/Melanoma_object.rds")

Receptors<-read.table("Receptor Gene Names.txt",sep="\t")

Receptors<-read.csv("/home/jtdecker/Receptor_interactions_curated.csv",header=TRUE)
Receptors<-Receptors[Receptors$Pair.Source=="known",]

library(mixOmics)
library(randomForest)

library(lars)
library(bnlearn)
library(minet)
library(igraph)
library(dplyr)
library(plyr)

pbmc<-TNBC

Receptor_names<-str_to_sentence(Receptors$Receptor.ApprovedSymbol)
Ligand_names<-str_to_sentence(Receptors$Ligand.ApprovedSymbol)

################ Linear Inference (Covariance) #########################################################

y<-viper_scores_df
x<-pbmc@assays$RNA@counts
x<-t(as.matrix(x[rownames(x)%in%Receptor_names,]))
x<-x[,colSums(x)>0]
LinearInference<-(cor((as.matrix(x)),y,method="pearson",use="pairwise.complete.obs"))####pearson correlation between values from t(n) and t(n+1)
write.table(LinearInference,"LI_Aaron_receptors.txt",sep="\t")

############### TSNI using LASSO constraints ############################################################
LASSO<-aaply(y,2,function(xx) {
  a<-lars(x,as.matrix(xx),type="lasso",use.Gram=TRUE) 
  m<-coef(a, s=.5, mode="lambda")
  m<- ifelse(m!=0,1,0)},.progress="text")
write.table(LASSO,"LASSO_TNBC_receptors.txt",sep="\t")

################# RANDOM FOREST ######################tm##################################################
RF<-aaply(y,2,function(xx) { 
  a<-randomForest(x,xx,mtry=round(sqrt(ncol(y))), ntree=100,  
                  importance=TRUE)
  m<-as.data.frame(importance(a))$IncNodePurity
  m<- ifelse(m>=quantile(m,0.9),1,0)},.progress="text")
write.table(RF,"RF_TNBC_receptors.txt",sep="\t")
################# Partial-Least Squares ################################################################
TDPLSR<-aaply(y,2,function(xx) { 
  a<-spls(x,xx,scale=FALSE,keepX=(round(.1*ncol(y))),ncomp =1)
  #m<-vip(a)
  m<- ifelse(a$loadings$X!=0,1,0)},.progress="text")
write.table(TDPLSR,"TDPLSR_TNBC_receptors.txt",sep="\t")

################# Mutual Information ###################################################################
table_for_analysis<-as.data.frame(cbind(x,y))
mim<-build.mim(table_for_analysis)
Mutual_Information<-aracne(mim)
Mutual_Information<-Mutual_Information[colnames(y),colnames(x)]
MI_mean<-aaply(as.matrix(colnames(y)),1,function(x) apply(Mutual_Information[which(rownames(Mutual_Information)==x),],2,sum))

############Network Assembly################################################


LI_conn<-read.table("LI_Aaron_receptors.txt")
yy<-LI_conn
LI_connections<-as.data.frame(ifelse(abs(yy)>=.3,1,0))
yy<-LI_connections
sum(LI_connections)

yy<-yy[rowSums(yy)>0,colSums(yy)>0]

yy<-yy[,colnames(yy)%in%(highly_variable_tfs)]


Network_Inference<-LI_connections

write.table(Network_Inference,"Aaron_scRNAseq_network.txt",sep="\t")

b<-colnames(Network_Inference)
a<-rownames(Network_Inference)

Edges<-matrix(0,ncol = 3)

for(tm in 1:length(b)){
  g<-b[tm]
  Edges<-rbind(Edges,aaply(as.matrix(a),1,function(x) c(x,g,Network_Inference[x,g])))}

Edges<-as.data.frame(Edges[Edges[,3]!="0",1:2])
colnames(Edges)<-c("From","To")
Nodes<-c(colnames(Network_Inference),rownames(Network_Inference))
Nodes<-Nodes[!duplicated(Nodes)]


net<-graph.data.frame(Edges,Nodes,directed=T)

net<-simplify(net,remove.multiple = F,remove.loops = T)
mm<- eigen_centrality(net)
quantile(mm$vector,0.9)
mm$vector
V(net)$type<-ifelse(mm$vector>quantile(mm$vector,0.9),2,1)
#colrs<-c(rgb(0,177/255,176/255))

#colrs<-c(rgb(0,177/255,176/255))
green<-c(rgb(0,177/255,176/255))
maize<-c(rgb(1,203/255,5/255))
blue=c("#00264c")
colrs<-c(green,maize)

V(net)$color<-colrs[V(net)$type]

deg<-(degree(net,mode="all"))



deg<-(degree(net))
V(net)$size<-sqrt(deg)
V(net)$color<-colrs[V(net)$type]
V(net)$frame.color<-colrs[V(net)$type]
V(net)$label.cex<-.1
V(net)$label.color<-"black"
V(net)$label.family<-"sans"
V(net)$label.font<-2



E(net)$width <- .01
E(net)$arrow.size<-.02
E(net)$color<-"black"


library(qgraph)
e <- get.edgelist(net,names=FALSE)


l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net),area=(vcount(net)^2),repulse.rad=(vcount(net)^1))
plot(net,layout=l,vertex.size=2)

########################Activity and Expression Assignments##########################################
# Network_Inference<-read.table("LI.txt")
# 
# melanoma_cell_types<-readRDS("/home/jtdecker/melanoma_cell_types.RDS")
# 
# Idents(pbmc) = as.factor(melanoma_cell_types)
# 
# TF_markers<-FindAllMarkers(pbmc,assay="dorothea")
# a<-colnames(Network_Inference)
# b<-rownames(Network_Inference)
# 
# Edges<-matrix(0,ncol = 3)
# 
# 
# Ligs<-pbmc@assays$RNA@scale.data
# 
# Ligs<-t(as.matrix(Ligs[rownames(Ligs)%in%Ligand_names,]))
# 
# 
# Rec_markers<-FindAllMarkers(pbmc,assay="RNA",features=colnames(x),only.pos=TRUE)
# TF_markers<-FindAllMarkers(pbmc,assay="dorothea",only.pos = TRUE)
# Lig_markers<-FindAllMarkers(pbmc,assay="RNA",features=colnames(Ligs),only.pos=TRUE)
# 
# 
# for(tm in 1:length(a)){
#   g<-a[tm]
#   Edges<-rbind(Edges,aaply(as.matrix(b),1,function(x) c(x,g,Network_Inference[x,g])))}
# 
# Edges<-as.data.frame(Edges[Edges[,3]!="0",1:2])
# colnames(Edges)<-c("From","To")
# 
# 
#   Network_score<-apply(as.matrix(unique(Edges$From)),1,function(xx){ ###Returns a list where each element is a transcription factor and a nxm matrix where n is a cell and m is a receptor
#    TF<-y[,colnames(y)==xx]
#   # TF<-TF-min(TF)
#    print(xx)
#     Rec<-x[,colnames(x)%in%Edges$To[Edges$From==xx] ]
#     Rec<-Rec-min(Rec)
#     if(length(Edges$To[Edges$From==xx]==1)){
#       Rec<-as.data.frame(Rec)
#       colnames(Rec)<-Edges$To[Edges$From==xx]}
#    pp<-TF*Rec
#     colnames(pp)<-paste(colnames(pp),"-",xx,sep="")
#     return(pp)
#     })
# 
# Network_score<-t(do.call(cbind.data.frame,Network_score))
# 
# Network_score<-apply(as.matrix(unique(Rec_markers$cluster)),1,function(xx){
#   TF_markers<-y[pbmc@active.ident=="T",]
#   
#   TFs<-TF_markers$gene[TF_markers$cluster==xx]
#   Recs<-Rec_markers$gene[Rec_markers$cluster==xx]
#   Ligs<-unique(Lig_markers$gene)
#   Connections<-Edges[Edges$From%in%TFs,]
#   Connections<-Connections[Connections$To%in%Recs,]
#   ww<-match(which(Receptors$Ligand.ApprovedSymbol%in%Ligs),which(Receptors$Receptor.ApprovedSymbol%in%Connections$To))
#   ww<-ww[!is.na(ww)]
#   Ligands<-which(Receptors$Receptor.ApprovedSymbol%in%Connections$To)[ww]
#   Receptor_Ligand<-data.frame(Receptor=Receptors$Receptor.ApprovedSymbol[Ligands],Ligand=Receptors$Ligand.ApprovedSymbol[Ligands])          
#   TF_Receptor<-Connections          
#   Ligands_CellTypes<-Lig_markers[Lig_markers$gene%in%Receptor_Ligand$Ligand,6:7]
#   
#   )
#   
#   
#   })
# 






names(Network_score)<-unique(Rec_markers$cluster)


cell_types<-read.csv("TNBC_scRNAseq_metadata.csv")


pbmc[["Network"]]<-CreateAssayObject(data=t(NN))



Idents(pbmc)<-cell_types$celltype_final
Network.markers<-FindAllMarkers(pbmc,assay="Network",min.pct = 0.3,only.pos = TRUE)

DefaultAssay(object = pbmc) <- "integrated"
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = rownames(pbmc), verbose = TRUE)
pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = TRUE)
pbmc <- FindClusters(pbmc,resolution=0.3, verbose = TRUE)

pbmc <- RunUMAP(pbmc, dims = 1:10, umap.method = "uwot", metric = "cosine")



pbmc<-NormalizeData(pbmc)
pbmc.markers <- FindAllMarkers(pbmc, only.pos=TRUE,logfc.threshold = 0.25, verbose = TRUE)
write.table(pbmc.markers,"TNBC_markers.txt",sep="\t")



DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5,group.by = "sample") + NoLegend()

Idents(pbmc)<-pbmc@meta.data$seurat_clusters


saveRDS(pbmc,"TNBC_scRNAseq_combined.rds")





cells<-cell_types[cell_types$NAME%in%rownames(x),c("NAME","celltype_final")]




x<-(pbmc@assays$RNA@counts)
x<-t(as.matrix(x[rownames(x)%in%Receptor_names,]))
x<-x[,colSums(x)>0]


#################Assemble Cell Type Networks######################

Network_score<-alply(as.matrix(c(0:15)),1,function(xx){
  print(xx)
  TF_markers<-y[pbmc$seurat_clusters==xx,]
  Rec_markers<-x[pbmc$seurat_clusters==xx,]

  sds<-apply(y,2,sd)

  TFs<-ifelse(TF_markers>(colMeans(y)+sds),1,0)
  Recs<-ifelse(Rec_markers>0,1,0)

  Edge_names<-paste(Edges$To,"-",Edges$From,sep="")

  NET<-alply(as.matrix(c(1:ncol(TFs))),1,function(xx){ ###Returns a list where each element is a transcription factor and a nxm matrix where n is a cell and m is a receptor
    TF<-TFs[,xx]
    pp<-TF*Recs
    #print(xx)
    colnames(pp)<-paste(colnames(TFs)[xx],"-",colnames(pp),sep="")
    return(pp)})


  NN<-(do.call(cbind,NET))
  NN<-NN[,colSums(NN)>0]
  NN<-NN[,colnames(NN)%in%Edge_names]
  return(NN)})

TF_Receptor<-Network_score

# 
# reg<-regulon[regulon$target%in%Ligand_names,] #defines which regulon exists for Ligand Names
# 
# Lig<-pbmc@assays$RNA@counts
# Lig<-t(as.matrix(Lig[rownames(Lig)%in%Ligand_names,]))
# 
# Network_score<-alply(as.matrix(c(0:15)),1,function(xx){
#   TF_markers<-y[pbmc$seurat_clusters==xx,]
#   Lig_markers<-Lig[pbmc$seurat_clusters==xx,]
#   
#   TFs<-ifelse(TF_markers>0,1,0)
#   Ligs<-ifelse(Lig_markers>0,1,0)
#   
#   Edge_names<-paste(reg$tf,"-",reg$target,sep="")
#   
#   NET<-alply(as.matrix(c(1:ncol(TFs))),1,function(xx){ ###Returns a list where each element is a transcription factor and a nxm matrix where n is a cell and m is a receptor
#     TF<-TFs[,xx]
#     pp<-TF*Ligs
#     #print(xx)
#     colnames(pp)<-paste(colnames(pp),"-",colnames(TFs)[xx],sep="")
#     return(pp)})
#   
#   
#   NN<-(do.call(cbind,NET))
#   NN<-NN[,colSums(NN)>0]
#   NN<-NN[,colnames(NN)%in%Edge_names]
#   return(NN)})
# 
# TF_Ligand<-Network_score
# 
# names(TF_Ligand)<-unique(cells$celltype_final)
# names(TF_Receptor)<-unique(cells$celltype_final)
# 
# Full_net<-cbind(NN,TF_Receptor)
# 


###Generate edges for each cluster for plotting###################### I REMOVED LIGANDS
Network_score_Edges<-alply(as.matrix(c(0:15)),1,function(xx){
  print(xx)
 # m<-data.frame(regulon$target,regulon$tf)
  #Excluded_connections<-paste(m[,1],m[,2],sep="-")
  
Tfh_R<-TF_Receptor[[xx+1]]
TR_sums<-colSums(Tfh_R)
Tfh_R<-Tfh_R[,TR_sums>0.1*nrow(Tfh_R)]
#Tfh_R<-Tfh_R[,!(colnames(Tfh_R)%in%Excluded_connections)]

#Tfh_L<-TF_Ligand[[xx+1]]
#TR_sums<-colSums(Tfh_L)
#Tfh_L<-Tfh_L[,TR_sums>0.1*nrow(Tfh_L)]

TFs<-data.frame(gene=unique(regulon$tf),type="TF")

#LL<-data.frame(gene=unique(Ligand_names),type="Ligand")
RR<-data.frame(gene=unique(Receptor_names),type="Receptor")
types<-rbind(TFs,RR,LL)
rownames(types)<-types[,1]

vv<-(strsplit(colnames(Tfh_R),"-"))
vv<-(do.call(rbind,vv))
vv<-vv[,1:2]
# vvv<-(strsplit(colnames(Tfh_L),"-"))
# vvv<-(do.call(rbind,vvv))
# vvv<-vvv[,1:2]
# vv<-rbind(vv,vvv)
#pp<-intersect(unique(vv[,1]),unique(vv[,2]))

#hh<-vv[vv[,1]%in%pp,]
#hhh<-vv[vv[,2]%in%pp,]
#vv<-rbind(hh,hhh)


Nodes<-c(vv[,1],vv[,2])
Nodes<-Nodes[!duplicated(Nodes)]

vv<-as.data.frame(vv)

return(vv)})






####Generate Nodes for the network for each cluster for plotting########################



Network_score_Nodes<-alply(as.matrix(c(0:15)),1,function(xx){
  m<-data.frame(regulon$target,regulon$tf)
  Excluded_connections<-paste(m[,1],m[,2],sep="-")
  Tfh_R<-TF_Receptor[[xx+1]]
  TR_sums<-colSums(Tfh_R)
  Tfh_R<-Tfh_R[,TR_sums>0.1*nrow(Tfh_R)]
  Tfh_R<-Tfh_R[,!(colnames(Tfh_R)%in%Excluded_connections)]
  
 # Tfh_L<-TF_Ligand[[xx+1]]
  #TR_sums<-colSums(Tfh_L)
  #Tfh_L<-Tfh_L[,TR_sums>0.1*nrow(Tfh_L)]
  
  
  TFs<-data.frame(gene=unique(regulon$tf),type="TF")
  
 # LL<-data.frame(gene=unique(Ligand_names),type="Ligand")
  RR<-data.frame(gene=unique(Receptor_names),type="Receptor")
  types<-rbind(TFs,RR)
  rownames(types)<-types[,1]
  
  vv<-(strsplit(colnames(Tfh_R),"-"))
  vv<-(do.call(rbind,vv))
  vv<-vv[,1:2]
  #vvv<-(strsplit(colnames(Tfh_L),"-"))
  #vvv<-(do.call(rbind,vvv))
  #vvv<-vvv[,1:2]
  #vv<-rbind(vv,vvv)
  #pp<-intersect(unique(vv[,1]),unique(vv[,2]))
  
  #hh<-vv[vv[,1]%in%pp,]
  #hhh<-vv[vv[,2]%in%pp,]
  #vv<-rbind(hh,hhh)
  
  
  Nodes<-c(vv[,1],vv[,2])
  Nodes<-Nodes[!duplicated(Nodes)]
  
  Nodes<-as.data.frame(Nodes,types)
  
  return(Nodes)})

########################################################################################

net<-graph.data.frame(Network_score_Edges[[7]],Network_score_Nodes[[7]],directed=T)


Nodes<-Network_score_Nodes[[7]]

TFs<-data.frame(gene=unique(regulon$tf),type="TF")

LL<-data.frame(gene=unique(Ligand_names),type="Ligand")
RR<-data.frame(gene=unique(Receptor_names),type="Receptor")
types<-rbind(TFs,LL,RR)
rownames(types)<-types[,1]


net<-simplify(net,remove.multiple = F,remove.loops = T)
mm<- eigen_centrality(net)
quantile(mm$vector,0.9)
mm$vector
V(net)$type<-types[match(Nodes[,1],types[,1]),2]
#colrs<-c(rgb(0,177/255,176/255))

#colrs<-c(rgb(0,177/255,176/255))
green<-c(rgb(0,177/255,176/255))
maize<-c(rgb(1,203/255,5/255))
blue=c("#00264c")
colrs<-c(green,maize,"red")

#V(net)$color<-"white"
#V(net)$label.color<-"white"
#E(net)$color<-"white"

deg<-(degree(net,mode="all"))



deg<-(degree(net))
V(net)$size<-sqrt(deg)
V(net)$color<-colrs[V(net)$type]
V(net)$frame.color<-colrs[V(net)$type]
V(net)$label.cex<-sqrt(deg)/20
V(net)$label.color<-"black"
V(net)$label.family<-"sans"
V(net)$label.font<-2



E(net)$width <- .01
E(net)$arrow.size<-.02
E(net)$color<-"black"


library(qgraph)
e <- get.edgelist(net,names=FALSE)


l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net),area=(vcount(net)^1.5),repulse.rad=(vcount(net)^2))
plot(net,layout=l)





######## Identify secreted ligands present in at least 50% of cells in the cluster###########
Secreted_Ligands<-alply(as.matrix(c(0:15)),1,function(xx){
  
  Lig_markers<-Lig[pbmc$seurat_clusters==xx,]
  
  Ligs<-ifelse(Lig_markers>0,1,0)
 
  Cell_number<-nrow(Lig_markers)
  Ligs<-colSums(Ligs)
  
  Secreted_Ligands<-names(Ligs[Ligs>.25*length(Ligs)])
  
  return(Secreted_Ligands)})


###Identify receptors in top 50% of nodes in the network receptor-tf-ligand network#####
library(igraph)
Important_Receptors<-alply(as.matrix(c(1:16)),1,function(xx){
  
  net<-graph.data.frame(Network_score_Edges[[xx]],Network_score_Nodes[[xx]],directed=T)
  
  
  Nodes<-Network_score_Nodes[[xx]]
  
  TFs<-data.frame(gene=unique(regulon$tf),type="TF")
  
  LL<-data.frame(gene=unique(Ligand_names),type="Ligand")
  RR<-data.frame(gene=unique(Receptor_names),type="Receptor")
  types<-rbind(TFs,LL,RR)
  rownames(types)<-types[,1]
  
  
  net<-simplify(net,remove.multiple = F,remove.loops = T)
  mm<- eigen_centrality(net)
  eigen<-mm$vector
  eigen<-eigen[names(eigen)%in%RR[,1]]
 # kk<-quantile(eigen,0.76)
  #eigen<-eigen[eigen>kk]
  
  return(eigen)})
  

######Find cell types that secrete factors from ealier threshold to important receptors in each cluster###########
#
#          We could put in the names of cell types here if we wanted to edit this code
#
#

library(rlist)

Receptor_ligand<-Receptors[,c(2,4)]

Receptor_ligand<-data.frame(str_to_sentence(Receptor_ligand$Ligand.ApprovedSymbol),str_to_sentence(Receptor_ligand$Receptor.ApprovedSymbol))

Receptor_Ligand_Connections<-alply(as.matrix(c(1:16)),1,function(xx){

  jj<-Important_Receptors[[xx]]

  RL<-Receptor_ligand[Receptor_ligand[,2]%in%names(jj),]

  qq<-list.search(Secreted_Ligands,.%in%RL[,1])

  NET<-alply(as.matrix(c(1:length(qq))),1,function(xx){ ###Returns a list where each element is a transcription factor and a nxm matrix where n is a cell and m is a receptor
    Possible_Ligs<-Secreted_Ligands[[xx]]
    Possible_Ligs<-Possible_Ligs[qq[[xx]]]
    Possible_connectins<-RL[RL[,1]%in%Possible_Ligs,]
    Edges<-paste(xx,"-",Possible_connectins[,1],"-",Possible_connectins[,2],sep="")
    return(Edges)})
  
  
  NN<-(do.call(cbind,NET))
  NN<-unique(unlist(as.data.frame(NN)))
  vv<-(strsplit(NN,"-"))
  vv<-(do.call(rbind,vv))
  vv<-vv[vv[,1]!=vv[,3],]
  vvv<-rbind(vv[,c(1:2)],vv[,c(2,3)])
  return(vvv)})

##########Add in TF-cluster marker connections#################
pbmc.markers<-read.table("TNBC_markers.txt",sep="\t")


Activated_markers<-alply(as.matrix(c(0:10)),1,function(xx){
  pp<-pbmc.markers[pbmc.markers$cluster==xx,]
  pp<-pp[pp$gene%in%regulon$target,c(2,7)]
  
  
  w<-regulon[regulon$target%in%pp[,2],]
  #w<-w[sign(w$mor)==sign(pp[1,1]),]
  TFs<-Network_score_Nodes[[xx+1]]
  w<-w[w$tf%in%t(TFs),]
  w<-w[w$mor>0,]
  ww<-data.frame(V1=w$tf,V2=w$target)
  return(data.frame(ww))})



Networks<-alply(as.matrix(c(0:15)),1,function(xx){
  w<-rbind(Network_score_Edges[[xx+1]],Receptor_Ligand_Connections[[xx+1]])#,Activated_markers[[xx+1]])
  return(w)})


############Plot Networks##############
Edge<-Networks[[8]]
Nodes<-unique(unlist(Networks[[8]]))

Edges<-as.matrix(Edge[,1])

Nodes[Nodes==1]<-"Monocytes"
Nodes[Nodes==2]<-"Macrophage"
Nodes[Nodes==3]<-"CD8 T cell"
Nodes[Nodes==4]<-"Monocyte"
Nodes[Nodes==5]<-"Fibroblast"
Nodes[Nodes==6]<-"Plasma"
Nodes[Nodes==7]<-"B cell"
Nodes[Nodes==8]<-"Breast Stromal"
Nodes[Nodes==9]<-"Cycling"
Nodes[Nodes==10]<-"Plasmacytoid DC"
Nodes[Nodes==11]<-"Mast"


Edges[Edges==1]<-"CD4+ T"
Edges[Edges==2]<-"CD8+ T"
Edges[Edges==3]<-"Tumor"
Edges[Edges==4]<-"Monocyte"
Edges[Edges==5]<-"Fibroblast"
Edges[Edges==6]<-"Plasma"
Edges[Edges==7]<-"B cell"
Edges[Edges==8]<-"Breast Stromal"
Edges[Edges==9]<-"Cycling"
Edges[Edges==10]<-"Plasmacytoid DC"
Edges[Edges==11]<-"Mast"

Edges<-data.frame(Edges,Edge[,2])

net<-graph.data.frame(Edges,Nodes,directed=T)
net<-simplify(net,remove.multiple = T,remove.loops = T)




TFs<-data.frame(gene=unique(regulon$tf),type="TF")

LL<-data.frame(gene=unique(Ligand_names),type="Ligand")
RR<-data.frame(gene=unique(Receptor_names),type="Receptor")


types<-rbind(TFs,LL,RR)
rownames(types)<-types[,1]


net<-simplify(net,remove.multiple = F,remove.loops = T)
mm<- eigen_centrality(net)
quantile(mm$vector,0.9)
mm$vector
V(net)$type<-3
V(net)$type<-types[match(Nodes,types[,1]),2]
V(net)$type[is.na(V(net)$type)]<-4
#colrs<-c(rgb(0,177/255,176/255))

#colrs<-c(rgb(0,177/255,176/255))
green<-c(rgb(0,177/255,176/255))
maize<-c(rgb(1,203/255,5/255))
blue=c("#00264c")
colrs<-c(green,maize,"red","lightblue")

#V(net)$color<-"white"
#V(net)$label.color<-"white"
#E(net)$color<-"white"

deg<-(degree(net,mode="all"))



deg<-(degree(net))
V(net)$size<-sqrt(deg)
V(net)$color<-colrs[as.factor(V(net)$type)]
V(net)$frame.color<-colrs[as.factor(V(net)$type)]
V(net)$label.cex<-sqrt(deg)/20
V(net)$label.color<-colrs[as.factor(V(net)$type)]
V(net)$label.family<-"sans"
V(net)$label.font<-2



E(net)$width <- .01
E(net)$arrow.size<-.02
E(net)$color<-"black"


library(qgraph)
e <- get.edgelist(net,names=FALSE)


l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(net),area=(vcount(net)^1.5),repulse.rad=(vcount(net)^2))
plot(net,layout=l)





new.cluster.ids<-c("CD4+ T","CD8+ T","Tumor","Monocyte","Fibroblast","Plasma","B cell","Breast Stroma","Cycling","Plasmacytoid DC","Mast cell")


names(new.cluster.ids)<-levels(pbmc)

pbmc<-RenameIdents(pbmc,new.cluster.ids)

UMAPPlot(pbmc,label=TRUE)



macrophages<-SplitObject(pbmc,split.by = "ident")

#####









































#
#
#






nodes_of_interest <- c("STAT4","TBX21")

# select the nodes having these namesV()
selnodes <- V(net)[name %in% nodes_of_interest]
# get their network neighborhood 
selegoV <- ego(net, order=1, nodes = selnodes, mode = "all", mindist = 0)

# turn the returned list of igraph.vs objects into a graph
G <- induced_subgraph(net,unlist(selegoV))

########################## plot the subgraph by greying out the stuff we don't care about##################


Edge_list<-as_edgelist(G)
ed<-unique(Edge_list[,1])
ed<-c(ed,unique(Edge_list[,2]))
in1<-incident_edges(net,nodes_of_interest,mode="all")

ed<-Nodes[!Nodes%in%ed]

V(net)[ed]$color<-"grey"
V(net)[ed]$frame.color<-"grey"
V(net)[ed]$label.color<-"grey"
E(net)[unlist(in1)]$color<-"black"

#E(net)[E(G)]$color<-"grey"
plot(net,layout=l)






















#










#
#



layout_new<-l[which(Nodes%in%ed),]
#l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(selegoG),area=(vcount(selegoG)^1.5),repulse.rad=(vcount(selegoG)^2))

plot(selegoG,layout=layout_new)



plot(selegoG,vertex.label=V(selegoG)$name)


All_edges<-cbind(Tfh_L,Tfh_R)

cell_type_edges<-paste(Edge_list[,1],"-",Edge_list[,2],sep="")

All_edges<-All_edges[,cell_type_edges]

cells_per<-rowSums(All_edges)

fre<-table(cells_per)

pie<-data.frame(Cell_Type = c("TH2","Other"),value = c(sum(fre)-fre[1],fre[1]))

ggplot(pie, aes(x="", y=value, fill=Cell_Type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_void() 




mmm<-table(pbmc$sample,pbmc$seurat_clusters)

####barplot function looks nicer

plot(mmm,color=c("red","orange","yellow","green","blue","purple"),main="Proportion of Cells in Each Cluster",cex=1.2)

clusters<-c("CD8 T","Plasma","Fibroblast","Macrophage complement","Immature T","Treg","Tumor1","B","Immature CD8 T","Cycling Tumor","Immature T_2","Macrohpage S100A8","Endothelial","Cycling T","Perivascular","Plasmacytoid DC","Keratincoyte")


mm<-pbmc
RenameIdents(mm,"0"="CD8 T","1"="Plasma","2"="Fibroblast","3"="Macrophage complement",
             "4"="Immature T","5"="Treg","6"="Tumor1","7"="B","8"="Immature CD8 T","9"="Cycling Tumor","10"="Immature T_2","11"="Macrohpage S100A8",
             "12"="Endothelial","13"="Cycling T","14"="Perivascular","15"="Plasmacytoid DC","16"="Keratincoyte")
Idents(mm)<-clusters
par(las=2)
plot(mmm,color=c("red","orange","yellow","green","blue","purple"),main="Proportion of Cells in Each Cluster",cex=1.2)


#######################################CellChat
library("CellChat")




data.input = pbmc@assays$RNA@counts # normalized data matrix
data.input<-normalizeData(data.input)
identity = data.frame(group = pbmc@active.ident, row.names = names(pbmc@active.ident)) # create a dataframe consisting of the cell labels
unique(identity$group) # check the cell labels
cellchat <- createCellChat(data = data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat<-addMeta(cellchat,meta=pbmc@meta.data)
cellchat <- setIdent(cellchat, ident.use = "seurat_clusters") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
cellchat@DB <- CellChatDB.use # set the used database in the object

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
#future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

saveRDS(cellchat,"TNBC cellchat.rds")






######################################  MONOCLE#######################################################

library(Seurat)
library(monocle3)

# The tutorial shown below and on subsequent pages uses two additional packages:
library(ggplot2)
library(dplyr)

pbmc<-readRDS("TNBC_scRNAseq_combined.rds")
pbmc<-subset(pbmc,ident=c("CD4+ T","CD8+ T"))


expression_matrix <- pbmc@assays$RNA@counts
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)






T_Cells<-subset(TNBC,idents=c("0","1"))

ifnb.list <- SplitObject(T_Cells, split.by = "patient")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:10,k.filter=5,k.score=5)
immune.combined <- IntegrateData(anchorset = immune.anchors,k.weight=10)


T_Cells<-immune.combined
DefaultAssay(T_Cells)<-"RNA"

T_Cells<-FindVariableFeatures(T_Cells)
T_Cells<-ScaleData(T_Cells)
T_Cells<-RunPCA(T_Cells, npcs = 30, verbose = TRUE)
# t-SNE and Clustering
T_Cells <- RunUMAP(T_Cells, reduction = "pca", dims = 1:20)
T_Cells <- FindNeighbors(T_Cells, reduction = "pca", dims = 1:20)
T_Cells<-FindClusters(T_Cells,resolution=.3)

UMAPPlot(T_Cells)


UMAPPlot(T_Cells,group.by="patient",label=TRUE)







#










####Visualization
pathways.show <- c("TGFb") 
vertex.receiver = seq(1,9) # a numeric vector
# Hierarchy plot
netVisual_aggregate(cellchat[cellchat@meta$tumor=="T",], signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
# Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",thresh=0.05)
netAnalysis_contribution(cellchat, signaling = pathways.show)

cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netVisual_signalingRole(cellchat, signaling = pathways.show)

nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,height=10)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")#,c("CD4+ T","CD8+ T","Tumor","Monocyte","Fibroblast","Plasma","B cell","Breast Stroma","Cycling","Plasmacytoid DC","Mast cell")nal")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional")
netVisual_embeddingZoomIn(cellchat, type = "functional")

cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural")
netVisual_embeddingZoomIn(cellchat, type = "structural")



























### T cell subset

Tcell<-subset(lung,CellType=="T Cells")


Tcell <- ScaleData(Tcell)
Tcell <- RunPCA(Tcell, features = rownames(Tcell), verbose = TRUE)
Tcell <- FindNeighbors(Tcell, dims = 1:10, verbose = TRUE)
Tcell <- FindClusters(Tcell,resolution=0.3, verbose = TRUE)

Tcell <- RunUMAP(Tcell, dims = 1:10, umap.method = "uwot", metric = "cosine")



Tcell<-NormalizeData(Tcell)
Tcell.markers <- FindAllMarkers(Tcell, only.pos=TRUE,logfc.threshold = 0.5 , verbose = TRUE)
write.table(Tcell.markers,"TNBC_markers.txt",sep="\t")





