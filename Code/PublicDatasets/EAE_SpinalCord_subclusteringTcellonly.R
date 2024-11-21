#Public Data
# Biological aging of CNS-resident cells alters the clinical course and
# immunopathology of autoimmune demyelinating disease
# doi: https://doi.org/10.1172/jci.insight.158153
# GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200901



setwd("/Users/lailarad/Documents/BI_EAE/Public scRNAseq Data/")
path = "/Users/lailarad/Documents/BI_EAE/Public scRNAseq Data/"



library(Seurat, lib.loc = .libPaths()[2])
packageVersion("Seurat")
library(Matrix)
packageVersion("Matrix")
library(SeuratObject, lib.loc = .libPaths()[2])
packageVersion("SeuratObject")
library(future, lib.loc = .libPaths()[2])
packageVersion("future.apply")
library(future.apply, lib.loc = .libPaths()[2])
library(ggprism)

library(dplyr)
library(plyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(purrr)
library(tibble)

file = list.files(path=path, pattern = ".gz" )
file

#Load Data----
#countsData1 <- read.csv(file = "GSE200901_RNAseq_counts.csv.gz", header = TRUE, row.names = 1)
#metadata1 <- read.csv(file = "GSE200901_RNAseq_metadata.csv.gz", header = TRUE, row.names = 1)

countsData <- read.csv(file = "GSE200901_counts.csv.gz", header = TRUE, row.names = 1)
metadata <- read.csv(file = "GSE200901_metadata.csv.gz", header = TRUE, row.names = 1)
rownames(metadata) = str_replace_all(rownames(metadata), "-", ".")
EAE_Sc <- CreateSeuratObject(counts = countsData, project = "EAE_SC", 
                             min.cells = 3, min.features = 200,
                             meta.data = metadata)


EAE_Sc[["percent.mt"]] = PercentageFeatureSet(EAE_Sc, pattern = "^mt-")


VlnPlot(EAE_Sc, features = c("nCount_RNA","nFeature_RNA", "percent.mt"), group.by = "age_day")
EAE_Sc$og_seurat_clusters = EAE_Sc$seurat_clusters



unique(EAE_Sc$age_day)

Idents(EAE_Sc)

#EAE_Sc = subset(EAE_Sc, idents = c("Young.SC.d6", "Young.SC.d10"))
unique(EAE_Sc$age_day)

#Integrate and Scale ----
ifnb.list <- SplitObject(EAE_Sc, split.by = "age_day")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20, k.filter = 25)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
DefaultAssay(immune.combined) <- "integrated"

s.genes <- str_to_sentence(cc.genes$s.genes)
g2m.genes <- str_to_sentence(cc.genes$g2m.genes)
immune.combined <- CellCycleScoring(immune.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



immune.combined <- ScaleData(immune.combined, vars.to.regress=c("S.Score","G2M.Score"), verbose = TRUE)
immune.combined<-FindVariableFeatures(immune.combined)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = TRUE)

# UMAP and Clustering ----
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined<-FindClusters(immune.combined,resolution=.4)

EAE_Sc <-immune.combined
#saveRDS(EAE_Sc,file="EAE_Sc_subcluster.rds")

#Post Clustering----
EAE_Sc = readRDS("EAE_Sc.rds")

DefaultAssay(EAE_Sc) = "RNA"

unique(EAE_Sc$orig.ident)

Idents(EAE_Sc) = EAE_Sc$orig.ident
EAE_Sc = subset(EAE_Sc, idents = c("Young-SC-d6", "Young-SC-d10"))

EAE_Sc[[]]

Idents(EAE_Sc) = EAE_Sc$integrated_snn_res.0.4

DimPlot(EAE_Sc, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 5, repel = T)

# ID cells----
DotPlot(EAE_Sc, features = c(
                            "Cd3e", "Cd4", "Cd8a","Ctla4", "Gata3","Foxp3","Rorc",
                            "Ncr1", "Klrb1c", "Klre1", 
                            "Cd79a", "Ms4a1", "Cd19","Ighm", 
                            "Cd300a", "Ly6c2", "Cd14", "Ccr2", 
                             "Apoe", "Ms4a7", "Ccr5", "Fcgr1", "Adgre1", "Cd68",
                            "Csf1", "Edn1", "Siglecf",
                            "Tmem119", "Trem2", "P2ry12","Aif1", "Csf1r",
                            "S100a8","S100a9", "Camp", "Retnlg", "Ly6g",
                            "Ifitm1", "Siglech", "Klk1", "Clec9a",
                             "Itgax", "Cd86", "Cd209a",
                            "Bst2", "Ly6a", "Nrp1",
                            "Cacnb3", "Fscn1", "Tmem150c", 
                            "Fcer1a", "Elane", "Kit","Sdc1", "Slamf7", "Jchain",
                            "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa"
), group.by = "integrated_snn_res.0.4", scale = F) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red")


DotPlot(EAE_Sc, features =c(
  "Ptprc", #CD45
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Foxp3", "Trac", # Tcells
  "Il7r","Cd5","Rorc","Gata3", "Maff","Il23r", #ILC
  "Mcpt1", "Fcer1a", "Kit", #Mast Cells
  "Stmn1", "Rrm1", #Helper NK/proliferating
  "Xcl1","Ncr1", "Klre1", "Klri1", "Gzma", "Gzmb",#NK
  "Col4a1","Flt1","Tie1", "Col3a1", "Col1a1", "Dcn", #Stromal
  # Cd11b. cd11c.  dc sign.  cd103             f4/80
  "Itgam", "Itgax", "Cd209a", "Itgae","Adgre4","Adgre1","Cd80", "Clec9a", 
  "Tmem119", "Trem2", "P2ry12","Siglech", "Hexb",#Microglia
  "Il4i1", "Flt3", "Ly75",
  "Cd14", "Cd16","Lyz1", #mon
  "Hp", "Mapk1",#chemoattractant
  "Mki67", "Top2a", # proliferation
  "C1qa", "C1qb","C1qc","C3", "C1ra", "C2", #complement
  "Cx3cr1", "Cd68","Ms4a7",#Mac
  
  "Mrc1", "Msr1", "Fcgr3",  "Ccr3", "Cxcl2", "Il1a", "Il1b","Ccr2", #Mac
  "S100a8", "S100a9", "Ly6g","Mmp8", "Mmp9", #Neut
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Tnfrsf13b", "Igha", "Slamf7", "Cd44", "Prdm1", "Xbp1", "Mcm5",  #Plasma 
  "Pecam1", "Cdh5", "Cxcl12", "Ivns1abp", "Itm2a", "Car4", "Ackr1", "Nr2f2", "Vwf", "Vcam1",#ECs
  "Pdgfrb","Vtn", "Higd1b", "Cspg4", "S1pr3", "Mcam", "Baiap3", "Ehd3",
  "Gfap", "S100b", "Clu", "Aqp4", "Agt", "Glul", "Slc1a3" 
)
, scale = F, assay = "RNA"
) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism(base_family = "Arial", base_size = 8) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1))


#```{r neutrophils, fig.height=12, fig.width=15, fig.asp=0.5 , fig.align='default'}
DotPlot(EAE_Sc, features = c("S100a9", "Mmp9", "Ly6g", "Cd177", "Ltf")) + RotatedAxis() + 
scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

#```{r monocytes, fig.height=12, fig.width=15, fig.asp=0.5 , fig.align='default'}
DotPlot(EAE_Sc, features = c("Ly6c2", "Ccr2", "F10", "Plac8", "Thbs1", "Cd14"), scale =F) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")



#```{r macrophages, fig.height=12, fig.width=15, fig.asp=0.5 , fig.align='default'}
DotPlot(EAE_Sc, features = c("Ms4a7", "Pf4", "Fapb4", "Gpnmb", "Gpr84")) + RotatedAxis() + 
scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

#```{r microglia, fig.height=12, fig.width=15, fig.asp=0.5 , fig.align='default'}
DotPlot(EAE_Sc, features = c("Gpr84", "Ptgs1", "P2ry12", "Gpr34", "Lag3", "Siglech", "Cst7"), scale = F) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")



#```{r dividing myeloid, fig.height=12, fig.width=15, fig.asp=0.5 , fig.align='default'}
DotPlot(EAE_Sc, features = c("Top2a", "Mki67", "Cdk1", "Ccnb2")) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")

FeaturePlot(EAE_Sc, 
            reduction = "umap", 
            features = c("Fcer1a", "Kit", "Elane"), 
            #sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)



get_conserved <- function(cluster){
  FindConservedMarkers(EAE_Sc,
                       ident.1 = cluster,
                       grouping.var = "age_day",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene")  %>%
    #left_join(y = unique(annotations[, c("gene_name", "description")]),
    #          by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}


Idents(EAE_Sc)
conserved_markers <- map_dfr(0:17, get_conserved)

# Extract top 10 markers per cluster
top20 <- conserved_markers %>% 
  mutate(avg_fc = (`Young-d6_avg_log2FC` +`Young-d10_avg_log2FC`) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 20, 
        wt = avg_fc)


plots <- VlnPlot(EAE_Sc, features = c(
  "Itgam", "Itgax", "Cd209a", "Itgae","Adgre4","Adgre1","Cd80", "Clec9a"
), 
split.by = "age_day",
pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

plots <- VlnPlot(EAE_Sc, features = c(
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Trac"
), 
split.by = "age_day",
pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

plots <- VlnPlot(EAE_Sc, features = c(
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", "Trac"
), 
split.by = "age_day",
pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 2)

# Cluster 0 
# Cluster 1 Chemoattractant Monocyte HP CD14 MHC-II
# Cluster 2 CD4+ T Cells
# Cluster 3 
# Cluster 4 Neutrophil "S100a9", "Camp", "Retnlg", "Ly6g",
# Cluster 5 
# Cluster 6 B cells CD19 Ms4a1/CD20
# Cluster 7 Microglia Tmem119, Trem2, P2ry12 
# Cluster 8 CD8+ T Cells
# Cluster 9 DCs Clec9a CD11c/Itgax CD209a/DCsign
# Cluster 10 CD4+ T Cells
# Cluster 11 
# Cluster 12 
# Cluster 13 NK Klrb1c Klre1 Ncr1
# Cluster 14 
# Cluster 15 Il4i1+ DCs Il4i1 Flt3 Ly75 Il1b MHC-II
# Cluster 16 B cells CD19 Ms4a1/CD20

# Cluster 0 Neutrophil Lrg1 Klf2 S100a9 or Eos Siglecf E
# Cluster 1 Mon "Ly6c2", "Ccr2", "F10", "Plac8", "Thbs1"
# Cluster 2 CD4 T cells CD3 CD4
# Cluster 3 Mon Cd14 Cd68
# Cluster 4 Mac "Ms4a7", "Pf4", "Fapb4", "Gpnmb" Adgre1
# Cluster 5 Microglia Tmem119, Trem2, P2ry12 
# Cluster 6 Neutrophil "S100a9", "Camp", "Retnlg", "Ly6g",
# Cluster 7 Mon "Ly6c2", "Ccr2", "F10", "Plac8", "Thbs1", "Cd14"
# Cluster 8 CD8 T cells CD3
# Cluster 9 B cells CD19 Ms4a1/CD20 CD79a
# Cluster 10 Microglia Tmem119, Trem2, P2ry12 
# Cluster 11 DCs Cd86, Siglech, Klk1, Clec9a
# Cluster 12 CD 4 T cells CD3 CD4
# Cluster 13 Mon "Ly6c2", "Ccr2", "F10", "Plac8", "Thbs1", "Cd14"
# Cluster 14 NK Klrb1c Klre1 Ncr1
# Cluster 15 B cells CD19 Ms4a1/CD20
# Cluster 16 CD4 T cells CD3 CD4
# Cluster 17 Mast Cells Cpa3 Fcer1a and Neutrophils Elane Ctsg MPO

table(Idents(EAE_Sc),EAE_Sc$day_group)

EAE_Sc = RenameIdents(EAE_Sc, 
                        "0" = "Neutrophils",
                        "1" = "Monocytes",
                        "2" = "CD4 T Cells",
                        "3" = "Monocytes",
                        "4" = "Macrophages",
                        "5" = "Microglia",
                        "6" = "Neutrophils",
                        "7" = "Monocytes",
                        "8" = "CD8 T Cells",
                        "9" = "B Cells",
                        "10" = "Microglia",
                        "11" = "DCs",
                        "12" = "CD4 T Cells",
                        "13" = "Monocytes",
                        "14" = "NKs",
                        "15" = "B Cells",
                        "16" = "CD4 T Cells",
                        "17" = "Neutrophils"
                        
)

cells.use <- WhichCells(EAE_Sc, expression=(`Cpa3` > 1 & Fcer1a >1), idents = c("Neutrophils"))
EAE_Sc <- SetIdent(EAE_Sc, cells = cells.use, value = 'Mast Cells')

cells.use <- WhichCells(EAE_Sc, idents = c("Mast Cells"))

EAE_Sc = subset(EAE_Sc, cells = cells.use, invert = TRUE)

EAE_Sc$CellIds = as.character(EAE_Sc$CellIds)

Idents(EAE_Sc) = "CellIds"

table(EAE_Sc$CellIds, EAE_Sc$age_day)


DimPlot(EAE_Sc, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 5, repel = T,
        split.by = "day_group")


getPalette = colorRampPalette(brewer.pal(9, "Set1"))
unique(EAE_Sc$CellIds)

DimPlot(EAE_Sc, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 9, repel = T, 
        group.by = "CellIds",
        cols = getPalette(10)) +
  theme_prism(base_size = 16, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 18)) +
  labs( title = "UMAP of Cell Types") 
ggsave(file = "UMAP_9clusters_EAE_Sc.pdf",
       units = "in",width = 11, height = 11)

ggsave(file = "UMAP_9clusters_EAE_Sc.png",
       units = "in",width = 11, height = 11, dpi = 400)



CellTypeTable = as.data.frame(proportions(table(EAE_Sc$CellIds,
                                                EAE_Sc$day_group), 
                                          margin = 2))
CellTypeTable$Freq = CellTypeTable$Freq*100
levels(CellTypeTable$Var1)


#display.brewer.all()



CellTypeTable$Var2 = factor(CellTypeTable$Var2, levels = c("d6", "d10"))
#png(file = "CellProportions_24clusters.png",units = "in",width = 14, height = 13, res = 400)
ggplot(CellTypeTable,aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism( base_size = 16)+
  labs(x = "Condition", y= "Percent of Cells") + 
  scale_fill_manual(values = getPalette(length(unique(EAE_Sc$CellIds))), name = "Cell Type") + 
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  )
ggsave(file = "CellProportions_10clusters_EAE_Sc.png",
       units = "in",width = 11, height = 11, dpi = 400)
dev.off()


#Cluster IDs
DotPlot(EAE_Sc, features =c(
  "Cd3e", "Cd4", "Cd8a", "Cd8b1", # Tcells
  "Lrg1", "Klf2", "S100a8", "S100a9", "Camp", "Retnlg", "Ly6g","Elane","Ctsg", #Neut
  "Fcer1a", "Kit","Cpa3", #Mast cell
  "Xcl1","Ncr1", "Klre1", "Gzma", "Gzmb",#NK
  "Ly6c2", "Ccr2", "F10", "Plac8", "Thbs1", "Cd14", "Cd68", #Mon
  "Ms4a7", "Pf4", "Gpnmb", "Adgre1",#Mac
  "Tmem119", "Trem2", "P2ry12", #Microglia
  "Itgax", "Clec9a","Siglech", "Klk1",  #DCs
  "H2-DMb1", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa", #MHC-II
  "Cd19",  "Ms4a1","Ighd","Ighm", "Cd22", #Bcells
  "Jchain", "Sdc1", "Igha"#Plasma 
)
, group.by = "CellIds") + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  theme_prism(base_family = "Arial", base_size = 20) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1))

ggsave(file = "DotplotMarkergenes_9clusters.png",
       units = "in",width = 17, height = 8, dpi = 400)

DefaultAssay(EAE_Sc)<-"RNA"

EAE_Sc$day_group <- factor(x = EAE_Sc$day_group, levels = c('d6', 'd10'))

VlnPlot(EAE_Sc, features = c("Ccr1", "Ccr2", "Ccr5", 
                           "Ccl2", "Ccl3", "Ccl5",
                           "Ccl6", "Ccl9"
                           ), 
        slot = "data", split.by = "day_group") + theme(legend.position = 'right')

VlnPlot(EAE_Sc, features = c("Ccr2", 
                             "Ccl2"
), 
slot = "data", split.by = "day_group") + theme(legend.position = 'right')


#T cell clustering ----
Idents(EAE_Sc)
Tcell_data = subset(EAE_Sc,idents = c("CD4 T Cells", "CD8 T Cells"))
Tcell_data[[]]
Idents(Tcell_data) = Tcell_data[["age_day"]]

ifnb.list1 <- SplitObject(Tcell_data, split.by = "age_day")

ifnb.list1 <- lapply(X = ifnb.list1, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

immune.anchors1 <- FindIntegrationAnchors(object.list = ifnb.list1, dims = 1:20, k.filter = 25)
immune.combined1 <- IntegrateData(anchorset = immune.anchors1, dims = 1:20)
DefaultAssay(immune.combined1) <- "integrated"

s.genes <- str_to_sentence(cc.genes$s.genes)
g2m.genes <- str_to_sentence(cc.genes$g2m.genes)
immune.combined1 <- CellCycleScoring(immune.combined1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



immune.combined1 <- ScaleData(immune.combined1, vars.to.regress=c("S.Score","G2M.Score", "percent.mt"), verbose = TRUE)
immune.combined1<-FindVariableFeatures(immune.combined1)
immune.combined1 <- RunPCA(immune.combined1, npcs = 30, verbose = TRUE)

# UMAP and Clustering 
immune.combined1 <- RunUMAP(immune.combined1, reduction = "pca", dims = 1:20)
immune.combined1 <- FindNeighbors(immune.combined1, reduction = "pca", dims = 1:20)
immune.combined1<-FindClusters(immune.combined1,resolution=.3)


Tcell_data <-immune.combined1


DimPlot(Tcell_data, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T)
DimPlot(Tcell_data, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 9, repel = T, group.by = "age_day")

Tcell.markers<-FindAllMarkers(Tcell_data,verbose = TRUE,only.pos = TRUE)

top20Tcell <- Tcell.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, 
        wt = avg_log2FC)

DefaultAssay(Tcell_data) = "RNA"

Tcellmarkers = c(
  "Cd3e","Cd3g", "Cd3d" ,"Cd8a","Cd8b1" ,"Cd4",
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
  "Il17a",  "Il17f" , "Il17re", "Il17rc", "Il17ra" ,"Il17rd", "Il17rb", "Il17d" ,
  "Klrb1c","Il6ra","Tgfb1","Il23a",
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
  "Trgv2", "Tcrg-C1", "Tcrg-C2","Tcrg-C4","Trdc", "Trdv4", "Trac", 
  "Klrg1", #CD8 memory
  "Klk1", "Nkg7", "Klrc2", "Ncr1", "Klre1", "Itgam", "Itgax",
  "Ftl1", "Fth1", "Apoe", "Cxcl10", "Cd274"
)

DotPlot(Tcell_data, features = Tcellmarkers, scale = F) + 
  RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "grey", high = "red") +  
  theme_prism(base_family = "Arial", base_size = 12) +
  theme(axis.text.x =element_text(angle = 45, hjust=1, vjust = 1))
dev.off()

table(Tcell_data$seurat_clusters, Tcell_data$age_day)

rownames(Tcell_data)
rownames(Tcell_data)[grepl("^Il7", rownames(Tcell_data))]

rownames(Tcell_data)[grepl("^Tcr", rownames(Tcell_data))]
rownames(Tcell_data)[grepl("^Trac", rownames(Tcell_data))]

rownames(Tcell_data)[grepl("^Tcrg", rownames(Tcell_data))]

VlnPlot(Tcell_data, features = c("Cd3e","Cd3g", "Cd3d","Cd8a","Cd4" ,"Foxp3", "Itgam", "Itgax",
                                 "Ftl1", "Fth1" ), split.by = "age_day")

VlnPlot(Tcell_data, features = c("Il17a","Il23r", "Rorc","Rora","Tcrg-V6", "Trgv2", "Trdc", "Trdv4"
                                 ), split.by = "age_day")


VlnPlot(Tcell_data, features = c("Cd8a","Cd8b1","Cd4" ,"Tcrg-V4", "Tcrg-V6", "Tcrg-C1", "Tcrg-C2", "Tcrg-V1", "Tcrg-C4", 
                                 "Tcrg-V5"), split.by = "age_day")

VlnPlot(Tcell_data, features = c("Lag3", "Izumo1r","Cd274"), split.by = "age_day")

VlnPlot(Tcell_data, features = c("Cd8a","Cd8b1","Cd4", "H2-DMb2","H2-Eb1", "H2-Aa", "H2-Ab1", "H2-DMa"),
        split.by = "age_day")

VlnPlot(Tcell_data, features = c("Cd40lg", "Icos", "Bcl6", "Il21r","Itgam", "Itgax",
                                 "Ccr6", "Cxcr3", "Cxcr5"), split.by = "age_day")

cells.use <- WhichCells(Tcell_data, expression=`H2-Eb1` > 0 |`Cd4` > 0 , idents = c("5"))

table(Tcell_data$integrated_snn_res.0.3, Tcell_data$age_day)

# Cluster 0 is CD4 effector memory Cd44+ CD62L- (Sell-)
# Cluster 1 is gamma/delta T Cd4-Cd8- "Il23r", "Rorc","Rora","Tcrg-V6", "Trgv2", "Trdc", "Trdv4",
# Cluster 2 is Treg CD4+ Foxp3+ Ctla4+
# Cluster 3 is CD11b+ CD4+ T
# Cluster 4 is Cytotoxic CD8+ T Nkg7+Gmza/b+Ccl5+
# Cluster 5 is CD11b+ CD11c+ CD4+ Tfh "Cd40lg", "Icos", "Bcl6", "Il21r",
# Cluster 6 are PDL1 CD4 T cells
# Cluster 7 are Th17 Rorc (Rorgamma) Il17+ Il23r

# Cluster 0 is CD4 effector memory Cd44+ CD62L- (Sell-)
# Cluster 1 is CD4 effector memory Cd44+ CD62L- (Sell-)
# Cluster 2 is CD4 effector memory Cd44+ CD62L- (Sell-)
# Cluster 3 is Cytotoxic CD8 T Gmza/b+Ccl5+
# Cluster 4 is Naive CD8 T Ccr7 
# Cluster 5 is CD4 CD11b
# Cluster 6 are CD4 CD11b
# Cluster 7 are gd CD4- Cd8- Il17a "Tcrg-C1", "Tcrg-C2","Tcrg-C4", "Trdc", "Trdv4",
# Cluster 8 are Treg CD4+ Foxp3+ Ctla4+

Tcell_data = RenameIdents(Tcell_data, 
                          "0" = "CD4 Tem",
                          "1" = "CD4 Tem",
                          "2" = "CD4 Tem",
                          "3" = "Cytotoxic CD8 T",
                          "4" = "Naive CD8 T",
                          "5" = "CD11b+ CD4+ T",
                          "6" = "CD11b+ CD4+ T",
                          "7" = "gd",
                          "8" = "Treg"
)

Tcell_data$TcellIDs = Idents(Tcell_data)
EAE_Sc[[]]
# Generate a new column called sub_cluster in the metadata
EAE_Sc$sub_cluster <- as.character(Idents(EAE_Sc))

# Change the information of cells containing sub-cluster information
EAE_Sc$sub_cluster[Cells(Tcell_data)] <- as.character(Idents(Tcell_data))
DimPlot(EAE_Sc, group.by = "sub_cluster", label = T)
getpalette(13)
color.use = c("orange","#232f75","#54b0e4","#b2df8a","#f781bf","#a65629","#4daf4a","#984ea3","#377eb8","#1c9e77" , "#e4211c", "#bc9dcc", 
             "#e3be00")

DimPlot(EAE_Sc, reduction = "umap", label = TRUE, pt.size = 0.1, label.size = 9, repel = T, 
        group.by = "sub_cluster",
        cols = c("orange","#232f75","#54b0e4","#b2df8a","#f781bf","#a65629","#4daf4a","#984ea3","#377eb8","#1c9e77" , "#e4211c", "#bc9dcc",  
                 "#e3be00")) +
  theme_prism(base_size = 16, base_fontface = "bold") + 
  theme(legend.text = element_text(size = 18)) +
  labs( title = "UMAP of Cell Types") 
ggsave(file = "UMAP_clustersTsubcluster_EAE_Sc_colorcorrect.pdf",
       units = "in",width = 11, height = 11)
dev.off()
ggsave(file = "UMAP_clustersTsubcluster_EAE_Sc.png",
       units = "in",width = 11, height = 11, dpi = 400)


CellTypeTable = as.data.frame(proportions(table(EAE_Sc$sub_cluster,
                                                EAE_Sc$day_group), 
                                          margin = 2))
CellTypeTable$Freq = CellTypeTable$Freq*100
levels(CellTypeTable$Var1)


#display.brewer.all()



CellTypeTable$Var2 = factor(CellTypeTable$Var2, levels = c("d6", "d10"))
#png(file = "CellProportions_24clusters.png",units = "in",width = 14, height = 13, res = 400)
ggplot(CellTypeTable,aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_col(width = 0.5, color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  theme_prism( base_size = 16)+
  labs(x = "Condition", y= "Percent of Cells") + 
  scale_fill_manual(values = getPalette(length(unique(EAE_Sc$sub_cluster))), name = "Cell Type") + 
  theme(legend.text = element_text(size = 16), 
        axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
        #plot.title = element_text(hjust = 0.5)
  )
ggsave(file = "CellProportions_Tsubclusters_EAE_Sc.png",
       units = "in",width = 11, height = 11, dpi = 400)
ggsave(file = "CellProportions_Tsubclusters_EAE_Sc.pdf",
       units = "in",width = 11, height = 11, dpi = 400)
ggsave(file = "CellProportions_Tsubclusters_EAE_Sc.svg",
       units = "in",width = 11, height = 11, dpi = 400)
dev.off()


#CellChat ----
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

library(Seurat,lib.loc = .libPaths()[2])
packageVersion("Seurat")
library(dplyr)
library(pheatmap)
library(tidyverse)
library(ggdendro)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(data.table)
library(GOplot)
library(clusterProfiler)
library(dittoSeq)
library(clusterProfiler)
library(RColorBrewer)
library(CellChat)
library(ComplexHeatmap)

setwd("/Users/lailarad/Documents/BI_EAE/Public scRNAseq Data/T cell subcluster/")

#robj_path = "/Users/lailarad/Documents/BI_EAE/Public scRNAseq Data/EAE_Sc.rds"
save_path = "/Users/lailarad/Documents/BI_EAE/Public scRNAseq Data/T cell subcluster/"

# Load the R object and save it as a loom file.
#data = readRDS(robj_path)
data = UpdateSeuratObject(object = EAE_Sc)
data[[]]
Idents(data) = "sub_cluster"
DimPlot(data)

data@meta.data$cells = data@active.ident

data@meta.data$short = data@meta.data$cells
levels(data@meta.data$short)
levels(Idents(data))
Idents(data) = data@meta.data$short



data = RenameIdents(data, 
                    "Neutrophils" = "Neut",
                    "Monocytes" = "Mon",
                    
                    "Macrophages" = "Mac",
                    "Microglia" = "Microglia",
                    
                    "B Cells" = "B Cells",
                    "DCs" = "DCs",
                    "NKs" = "NKs"
                   
                    )
data@meta.data$short = data@active.ident
DimPlot(data, reduction = 'umap', label = T)

levels(data@meta.data$short)

Idents(data) = data@meta.data$day_group
D6 = subset(data, idents = 'd6')
D10 = subset(data, idents = 'd10')

#Day 6 cellchat analysis ----
Idents(D6) = D6@meta.data$short
D6 = SetIdent(D6, value = D6@meta.data$short)
levels(Idents(D6))
cellchat <- createCellChat(object = D6, group.by = 'short', assay = "RNA")
cell.labels = levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

# Set Receptor Ligand Database --------------------------------------------
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data

showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF…
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T…
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB…
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF…
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T…
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago…
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""…
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep…
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435…
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre…
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2…

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

view(CellChatDB$interaction)

# Preprocessing Expression Data for Cell-Cell Communication Analysis --------
# subset the expression data of signaling genes for saving computation cost
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

number = which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") # 1887
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-number,]
number = which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") #1900
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-number,]

# set the used database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.mouse, raw.use = F)

# Compute Communicaiton probability and infer network ---------------------
cellchat <- computeCommunProb(cellchat) #remove pop.size = F
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 4)
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-08-22 22:38:49.873196]"
# |=================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-08-22 22:41:50.0529]"
# Warning message:
#   In UseMethod("depth") :
#   no applicable method for 'depth' applied to an object of class "NULL"


# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# 

df.net$interaction_ST = paste0(df.net$interaction_name, "_",
                                        df.net$source, "_",
                                        df.net$target)

df.net.path <- subsetCommunication(cellchat, slot.name = "netP") # look at level of signaling pathways

df.net.ccl.cxcl <- subsetCommunication(cellchat, signaling = c("CCL", "CXCL"))

df.net.ccl.cxcl$interaction_ST = paste0(df.net.ccl.cxcl$interaction_name, "_",
                                        df.net.ccl.cxcl$source, "_",
                                        df.net.ccl.cxcl$target)


# Infer cell-cell communication at signaling pathway level ----------------
cellchat <- computeCommunProbPathway(cellchat)


# Calculate aggregated communication network ------------------------------

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(3,5), mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


cellchat@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)


netVisual_aggregate(cellchat, signaling = c("CCL"), layout = "chord")
netVisual_aggregate(cellchat, signaling = c("CXCL"), layout = "chord")
netVisual_bubble(cellchat,  signaling = c("CCL","CXCL"), remove.isolate = T)
netVisual_chord_gene(cellchat, signaling = c("CCL","CXCL"),legend.pos.x = 8, thresh = 0.01)

#png(file = "CCLCXCL_geneChord_MacMicrogliaSources_Day6.png",units = "in",width = 5.5, height = 5, res = 400)
netVisual_chord_gene(cellchat, signaling = c("CCL","CXCL"),legend.pos.x = 2, sources.use = c(4,3),
                     small.gap = 4)
dev.off()
#ggsave(file = "CCLCXCL_geneChord_MacMicrogliaSources.png", width = 4 , height = 4, dpi = 400)

png(file = "CCLCXCL_geneChord_all_Day6.png",units = "in",width = 5.5, height = 5, res = 400)
netVisual_chord_gene(cellchat, signaling = c("CCL","CXCL"),legend.pos.x = 2, #sources.use = c(4,5),
                     small.gap = 3, big.gap = 10)
dev.off()

netVisual_bubble(cellchat,  signaling = c("CCL","CXCL"), remove.isolate = T)
ggsave(file = "CCLCXCL_bubble_all_Day6.png", width = 16 , height = 6, dpi = 400)


netVisual_bubble(cellchat,  remove.isolate = T)
ggsave(file = "=allsignaling_bubble_all_Day6.png", width = 16 , height = 16, dpi = 400)

dev.off()
netAnalysis_contribution(cellchat, signaling = c("CCL","CXCL"), targets.use = c(4,5))


# Part IV: Systems analysis of cell-cell communication network ----
# Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
# (A) Compute and visualize the network centrality scores

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = "CCL", width = 8, height = 2.5, font.size = 10)

# (B) Visualize dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


saveRDS(cellchat, file = "cellchat_SC_Day6_subcluster.rds")
#Day 10 cellchat analysis ----
Idents(D10) = D10@meta.data$short
D10 = SetIdent(D10, value = D10@meta.data$short)
cellchat10 <- createCellChat(object = D10, group.by = 'short', assay = "RNA")
cell.labels = levels(cellchat10@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat10@idents)) # number of cells in each cell group
groupSize

# Set Receptor Ligand Database --------------------------------------------
# CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data

#showDatabaseCategory(CellChatDB)

# Show the structure of the database
#dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF…
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T…
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB…
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF…
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T…
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago…
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""…
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep…
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435…
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre…
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2…

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB



# Preprocessing Expression Data for Cell-Cell Communication Analysis --------
# subset the expression data of signaling genes for saving computation cost
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# number = which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") # 1887
# CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-number,]
# number = which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") #1900
# CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-number,]

# set the used database in the object
cellchat10@DB <- CellChatDB.use

cellchat10 <- subsetData(cellchat10) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat10 <- identifyOverExpressedGenes(cellchat10)
cellchat10 <- identifyOverExpressedInteractions(cellchat10)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.mouse, raw.use = F)

# Compute Communicaiton probability and infer network ---------------------
cellchat10 <- computeCommunProb(cellchat10) #remove pop.size = F
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat10 <- filterCommunication(cellchat10, min.cells = 4)
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-08-22 22:38:49.873196]"
# |=================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-08-22 22:41:50.0529]"
# Warning message:
#   In UseMethod("depth") :
#   no applicable method for 'depth' applied to an object of class "NULL"


# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
df10.net <- subsetCommunication(cellchat10) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

df10.net$interaction_ST = paste0(df10.net$interaction_name, "_",
                               df10.net$source, "_",
                               df10.net$target)

df10.net.path <- subsetCommunication(cellchat10, slot.name = "netP") # look at level of signaling pathways

df10.net.ccl.cxcl <- subsetCommunication(cellchat10, signaling = c("CCL", "CXCL"))

df10.net.ccl.cxcl$interaction_ST = paste0(df10.net.ccl.cxcl$interaction_name, "_",
                                          df10.net.ccl.cxcl$source, "_",
                                          df10.net.ccl.cxcl$target)

# Infer cell-cell communication at signaling pathway level ----------------
cellchat10 <- computeCommunProbPathway(cellchat10)


# Calculate aggregated communication network ------------------------------

cellchat10 <- aggregateNet(cellchat10)

groupSize <- as.numeric(table(cellchat10@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat10@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat10@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat10@net$weight
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(3,5), mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


cellchat10@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat10@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat10@idents)


netVisual_aggregate(cellchat10, signaling = c("CCL"), layout = "chord")
netVisual_aggregate(cellchat10, signaling = c("CXCL"), layout = "chord")
netVisual_bubble(cellchat10,  signaling = c("CCL","CXCL"), remove.isolate = T)
netVisual_chord_gene(cellchat10, signaling = c("CCL","CXCL"),legend.pos.x = 8, thresh = 0.01)
plotGeneExpression(cellchat10, signaling = c("CCL","CXCL"), enriched.only = TRUE, type = "violin")

par(mfrow=c(1,1))
netVisual_heatmap(cellchat10, signaling = c("CCL"), color.heatmap = "Reds")

table(cellchat10@idents)
netAnalysis_contribution(cellchat10, signaling = c("CCL","CXCL"))
netAnalysis_contribution(cellchat10, signaling = c("CCL","CXCL"), sources.use = c(2,7)) #Mon, Dcs
netAnalysis_contribution(cellchat10, signaling = c("CCL","CXCL"), targets.use = c(2,7))

netAnalysis_contribution(cellchat10, signaling = c("CCL","CXCL"), sources.use = c(4,5)) #Mac, Microglia
netVisual_chord_gene(cellchat10, signaling = c("CCL","CXCL"),legend.pos.x = 8, sources.use = c(4,5))
dev.off()
netAnalysis_contribution(cellchat10, signaling = c("CCL","CXCL"), targets.use = c(4,5))

pairLR.CCL.CXCL <- extractEnrichedLR(cellchat10, signaling = c("CCL","CXCL"), geneLR.return = FALSE)
LR.show = pairLR.CCL.CXCL[6,]
netVisual_individual(cellchat10, signaling = c("CCL","CXCL"), pairLR.use = LR.show, layout = "chord")

plotGeneExpression(cellchat10, signaling = c("CCL","CXCL"), enriched.only = TRUE, type = "violin")
ggsave(file = "CCLCXCLexpressionDay10.png", width = 16, height = 32, pointsize = 12)

netAnalysis_contribution(cellchat10, signaling = c("CCL","CXCL"), sources.use = c(4,5),
                         title = "Contribution of each L-R pair at Day 10 \n Macrophages and Microglia as sources") #Mac, Microglia
ggsave(file = "CCLCXCL_LRcontribution_MacMicrogliaSources_Day10.png", width = 4 , height = 4, dpi = 400)

dev.off()

png(file = "CCLCXCL_geneChord_MacMicrogliaSources_Day10.png",units = "in",width = 5.5, height = 5, res = 400)
netVisual_chord_gene(cellchat10, signaling = c("CCL","CXCL"),legend.pos.x = 2, sources.use = c(4,5),
                     small.gap = 4 )
dev.off()
#ggsave(file = "CCLCXCL_geneChord_MacMicrogliaSources.png", width = 4 , height = 4, dpi = 400)

png(file = "CCLCXCL_geneChord_all_Day10.png",units = "in",width = 5.5, height = 5, res = 400)
netVisual_chord_gene(cellchat10, signaling = c("CCL","CXCL"),legend.pos.x = 2, #sources.use = c(4,5),
                     small.gap = 3, big.gap = 10)
dev.off()

netVisual_bubble(cellchat10,  signaling = c("CCL","CXCL"), remove.isolate = T)
ggsave(file = "CCLCXCL_bubble_all_Day10.png", width = 16 , height = 6, dpi = 400)


netVisual_bubble(cellchat10,  remove.isolate = T)
ggsave(file = "allsignaling_bubble_all_Day10.png", width = 16 , height = 16, dpi = 400)

# Part IV: Systems analysis of cell-cell communication network ----
# Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
# (A) Compute and visualize the network centrality scores

# Compute the network centrality scores
cellchat10 <- netAnalysis_computeCentrality(cellchat10, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat10, signaling = "CCL", width = 8, height = 2.5, font.size = 10)

# (B) Visualize dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat10)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat10, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


saveRDS(cellchat10, file = "cellchat_SC_Day10_subcluster.rds")
#Compare Day 6 and Day 10 ----

object.list <- list(Day6 = cellchat, Day10 = cellchat10)
cellchatScDays <- mergeCellChat(object.list, add.names = names(object.list))

unique(cellchatScDays@meta$day_group)

gg1 <- compareInteractions(cellchatScDays, show.legend = F, group = c(1,2) )
gg2 <- compareInteractions(cellchatScDays, show.legend = F, group = c(1,2), measure = "weight",  )
gg1 + gg2

# two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased
# (or decreased) signaling in the second dataset compared to the first one.
#red more in Day 10
# blue less in Day 10

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchatScDays, weight.scale = T)
netVisual_diffInteraction(cellchatScDays, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchatScDays)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchatScDays, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, 
                   label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- rankNet(cellchatScDays, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchatScDays, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

object.list


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Day10"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchatScDays <- identifyOverExpressedGenes(cellchatScDays, group.dataset = "datasets", 
                                           pos.dataset = pos.dataset, features.name = features.name, 
                                           only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchatScDays, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in Day10
net.up <- subsetCommunication(cellchatScDays, net = net, datasets = "Day10",ligand.logFC = 0.2)
#net.down <- subsetCommunication(cellchatScDays, net = net, datasets = "Day6",ligand.logFC = 0.2)


#write.csv(net.up,paste(plotPath,"/net.up_Day7_MogvsOVA.csv", sep = ""), row.names = T)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in Day6, i.e.,downregulated in Day10
net.down <- subsetCommunication(cellchatScDays, net = net, datasets = "Day6",ligand.logFC = -0.2, receptor.logFC = -0.2)

#write.csv(net.down,paste(plotPath,"/net.down_Day7_MogvsOVA.csv", sep = ""), row.names = T)

gene.up <- extractGeneSubsetFromPair(net.up, cellchatScDays)
gene.down <- extractGeneSubsetFromPair(net.down, cellchatScDays)

pairLR.use.up = net.up[, "interaction_name", drop = F]
c(3,8,9,10,15)
gg1 <- netVisual_bubble(cellchatScDays, pairLR.use = pairLR.use.up, #targets.use = c(Tcells, Bcells), 
                        #sources.use = c(Tcells, Bcells), 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg1
ggsave(file = "UprrgulatedDay10comparedtoDay6.png", width = 16 , height = 16, dpi = 400)

pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchatScDays, pairLR.use = pairLR.use.down, #targets.use = c(Tcells, Bcells), 
                        #sources.use = c(Tcells, Bcells), 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


gg2
ggsave(file = "DownrgulatedDay10comparedtoDay6.png", width = 16 , height = 16, dpi = 400)
dev.off()
gg2 = netVisual_chord_gene(object.list[[2]], #sources.use = c(3,8,9,10,15),targets.use = c(3,8,9,10,15), 
                            big.gap = 1 ,slot.name = 'net', net = net.down, lab.cex = 0.5, 
                            small.gap = 1, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))



pathways.show <- c("CCL") 
#png(file = "CCL_Genechorddiagram_Upreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated in" ,names(object.list)[2]),
                     lab.cex = 1)
dev.off()

pairLR.use.cxcl9 = net.down[which(net.down$interaction_name == "CXCL9_CXCR3"),]
pathways.show <- c("CXCL") 
#png(file = "CCL_Genechorddiagram_Upreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.down, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in EAE spinal cord over time"),
                     lab.cex = 1)

pathways.show <- c("OSM") 
#png(file = "CCL_Genechorddiagram_Upreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.down, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in EAE spinal cord over time"),
                     lab.cex = 1)

#png(file = "CCL_Genechorddiagram_Downreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.down, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list)[2]),
                     lab.cex = 1)
dev.off()

pathways.show <- c( "LAMININ")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i], unique(cellchat@meta$Day)))
}


pathways.show <- c( "ITGAL-ITGB2") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i], unique(cellchat@meta$Day)))
}

dev.off()

pathways.show <- c( "ICAM") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list)[i], unique(cellchat@meta$Day)))
}

#net.up.scaf.Day9v7 = read.csv("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/D7vD9 Mog/net.up.Day9v7.csv")
#net.down.scaf.Day9v7 = read.csv("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/D7vD9 Mog/net.down.Day9v7.csv")


pathways.show <- c( "CCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list)[i], unique(cellchat@meta$Day)))
}

unique(common.net.up$pathway_name)
pathways.show <- c( "CCL") 
getwd()
png(file = "CCL_Genechorddiagram_Up_SC_overtime.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "CCL_Genechorddiagram_Up_SC_overtime.pdf",width = 8, height = 7)

netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated \nin the spinal cord over time" ),
                     lab.cex = 1)
dev.off()
saveRDS(cellchatScDays, file = "cellchat_SC_comp6_10.rds")

#Load All scaffold data ----
#setwd("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/")

robj_path = "/Users/lailarad/Documents/ProgEAE_scRNAseq/LMR_16renamedclusters_wTcellsubcluster.rds"
#save_path = "/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/"

# Load the R object and save it as a loom file.
data = readRDS(robj_path)
data = UpdateSeuratObject(object = data)
data[[]]
Idents(data) = "sub_cluster"
DimPlot(data)

data@meta.data$cells = data@active.ident

data@meta.data$short = data@meta.data$cells
levels(data@meta.data$short)
levels(Idents(data))
Idents(data) = data@meta.data$short

cellchat@idents
data = RenameIdents(data, 
                    'Monocytes' = 'Mon',
                    "Macrophage" = "Mac",
                    "Complement Macrophage" = "Mac",
                    "Neutrophil" = "Neut",
                    "NK Cell" = "NKs",
                    "Stromal Cell" = "Stromal",
                    "CD11b+ DC" = "DCs",
                    "Helper NK" = "NKs",
                    "B Cell" = "B Cells",
                    "CD103+ DC" = "DCs",
                    "Chemoattractant Monocyte" = "Mon",
                    "Il4i1+ DC" = "DCs",
                    "Plasmablast" = "Plasmablast",
                    "Naive CD4" = "Naive CD4 T",
                    "CD4+ Effector Memory" = "CD4 Tem",
                    "Naive CD8" = "Naive CD8 T",
                    "Cytotoxic CD8 Tem" = "Cytotoxic CD8 T",
                    "gd" = "gd",
                    "Treg" = "Treg",
                    "Ferritin Myeloid Cells" = "Fer Myeloid")

data@meta.data$short = data@active.ident
DimPlot(data, reduction = 'umap', label = T)

Idents(data) = data@meta.data$Day
D7 = subset(data, idents = 'D7')


Idents(D7) = D7@meta.data$Antigen
mog = subset(D7, idents = "MOG")
ova = subset(D7, idents = "OVA")

Idents(mog) = mog@meta.data$short
mog = SetIdent(mog, value = mog@meta.data$short)

#MOG Day 7 ----
cellchat.D7MOG <- createCellChat(object = mog, group.by = 'short', assay = "RNA")
cell.labels = levels(cellchat.D7MOG@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
groupSize

# set the used database in the object
cellchat.D7MOG@DB <- CellChatDB.use

cellchat.D7MOG <- subsetData(cellchat.D7MOG) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat.D7MOG <- identifyOverExpressedGenes(cellchat.D7MOG)
cellchat.D7MOG <- identifyOverExpressedInteractions(cellchat.D7MOG)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.mouse, raw.use = F)

# Compute Communicaiton probability and infer network ---------------------
cellchat.D7MOG <- computeCommunProb(cellchat.D7MOG) #remove pop.size = F
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.D7MOG <- filterCommunication(cellchat.D7MOG, min.cells = 4)
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-08-22 22:38:49.873196]"
# |=================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-08-22 22:41:50.0529]"
# Warning message:
#   In UseMethod("depth") :
#   no applicable method for 'depth' applied to an object of class "NULL"


# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
cellchat.D7MOG.net <- subsetCommunication(cellchat.D7MOG) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways


cellchat.D7MOG.net$interaction_ST = paste0(cellchat.D7MOG.net$interaction_name, "_",
                                                    cellchat.D7MOG.net$source, "_",
                                                    cellchat.D7MOG.net$target)

cellchat.D7MOG.net.path <- subsetCommunication(cellchat.D7MOG, slot.name = "netP") # look at level of signaling pathways

cellchat.D7MOG.net.ccl.cxcl <- subsetCommunication(cellchat.D7MOG, signaling = c("CCL", "CXCL"))
cellchat.D7MOG.net.ccl.cxcl$interaction_ST = paste0(cellchat.D7MOG.net.ccl.cxcl$interaction_name, "_",
                                                    cellchat.D7MOG.net.ccl.cxcl$source, "_",
                                                    cellchat.D7MOG.net.ccl.cxcl$target)


common_interactionsEarly = intersect(df.net.ccl.cxcl$interaction_ST, cellchat.D7MOG.net.ccl.cxcl$interaction_ST)

# Infer cell-cell communication at signaling pathway level ----------------
cellchat.D7MOG <- computeCommunProbPathway(cellchat.D7MOG)


# Calculate aggregated communication network ------------------------------

cellchat.D7MOG <- aggregateNet(cellchat.D7MOG)

groupSize <- as.numeric(table(cellchat.D7MOG@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.D7MOG@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.D7MOG@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat.D7MOG@net$weight
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(3,5), mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


cellchat.D7MOG@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat.D7MOG@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat.D7MOG@idents)

netVisual_chord_gene(cellchat.D7MOG, signaling = c("CCL","CXCL"),legend.pos.x = 2, #sources.use = c(4,5),
                     small.gap = 3, big.gap = 10)
dev.off()

# Part IV: Systems analysis of cell-cell communication network ----
# Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
# (A) Compute and visualize the network centrality scores

# Compute the network centrality scores
cellchat.D7MOG <- netAnalysis_computeCentrality(cellchat.D7MOG, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.D7MOG, signaling = "CCL", width = 8, height = 2.5, font.size = 10)

# (B) Visualize dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat.D7MOG)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat.D7MOG, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#OVA D7 ----
Idents(ova) = ova@meta.data$short
ova = SetIdent(ova, value = ova@meta.data$short)
ova[[]]
cellchat.D7OVA <- createCellChat(object = ova, group.by = 'short', assay = "RNA")
levels(cellchat.D7OVA@idents) # show factor levels of the cell labels
as.numeric(table(cellchat.D7OVA@idents)) # number of cells in each cell group


# set the used database in the object
cellchat.D7OVA@DB <- CellChatDB.use

cellchat.D7OVA <- subsetData(cellchat.D7OVA) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat.D7OVA <- identifyOverExpressedGenes(cellchat.D7OVA)
cellchat.D7OVA <- identifyOverExpressedInteractions(cellchat.D7OVA)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.mouse, raw.use = F)

# Compute Communicaiton probability and infer network ---------------------
cellchat.D7OVA <- computeCommunProb(cellchat.D7OVA) #remove pop.size = F
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.D7OVA <- filterCommunication(cellchat.D7OVA, min.cells = 4)
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-08-22 22:38:49.873196]"
# |=================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-08-22 22:41:50.0529]"
# Warning message:
#   In UseMethod("depth") :
#   no applicable method for 'depth' applied to an object of class "NULL"


# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
cellchat.D7OVA.net <- subsetCommunication(cellchat.D7OVA) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways


cellchat.D7OVA.net$interaction_ST = paste0(cellchat.D7OVA.net$interaction_name, "_",
                                                    cellchat.D7OVA.net$source, "_",
                                                    cellchat.D7OVA.net$target)

cellchat.D7OVA.net.path <- subsetCommunication(cellchat.D7OVA, slot.name = "netP") # look at level of signaling pathways

cellchat.D7OVA.net.ccl.cxcl <- subsetCommunication(cellchat.D7OVA, signaling = c("CCL", "CXCL"))
cellchat.D7OVA.net.ccl.cxcl$interaction_ST = paste0(cellchat.D7OVA.net.ccl.cxcl$interaction_name, "_",
                                                    cellchat.D7OVA.net.ccl.cxcl$source, "_",
                                                    cellchat.D7OVA.net.ccl.cxcl$target)

intersect(common_interactionsEarly, cellchat.D7OVA.net.ccl.cxcl$interaction_ST)

#common_interactionsEarly = intersect(df.net.ccl.cxcl$interaction_ST, cellchat.D7MOG.net.ccl.cxcl$interaction_ST)

# Infer cell-cell communication at signaling pathway level ----------------
cellchat.D7OVA <- computeCommunProbPathway(cellchat.D7OVA)


# Calculate aggregated communication network ------------------------------

cellchat.D7OVA <- aggregateNet(cellchat.D7OVA)

groupSize <- as.numeric(table(cellchat.D7OVA@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.D7OVA@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.D7OVA@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat.D7OVA@net$weight
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(3,5), mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


cellchat.D7OVA@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat.D7OVA@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat.D7OVA@idents)

netVisual_chord_gene(cellchat.D7OVA, signaling = c("CCL","CXCL"),legend.pos.x = 2, #sources.use = c(4,5),
                     small.gap = 3, big.gap = 10)
dev.off()

# Part IV: Systems analysis of cell-cell communication network ----
# Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
# (A) Compute and visualize the network centrality scores

# Compute the network centrality scores
cellchat.D7OVA <- netAnalysis_computeCentrality(cellchat.D7OVA, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.D7OVA, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# (B) Visualize dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat.D7OVA)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat.D7OVA, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#Compare Day 7 MOG v OVA ----

object.list3 <- list(MOG = cellchat.D7MOG, OVA = cellchat.D7OVA)
cellchatMvO7 <- mergeCellChat(object.list3, add.names = names(object.list3))

unique(cellchatMvO7@meta)

gg1 <- compareInteractions(cellchatMvO7, show.legend = F, group = c(1,2) )
gg2 <- compareInteractions(cellchatMvO7, show.legend = F, group = c(1,2), measure = "weight" )
gg1 + gg2

# two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased
# (or decreased) signaling in the second dataset compared to the first one.
#red more in Day 10
# blue less in Day 10

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchatMvO7, weight.scale = T)
netVisual_diffInteraction(cellchatMvO7, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchatMvO7)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchatMvO7, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list3, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list3)) {
  netVisual_circle(object.list3[[i]]@net$count, weight.scale = T, 
                   label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list3)[i]))
}


num.link <- sapply(object.list3, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list3)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list3[[i]], 
                                               title = names(object.list3)[i],
                                               weight.MinMax = weight.MinMax) + xlim(0,12) + ylim(0,10)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- rankNet(cellchatMvO7, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchatMvO7, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

object.list3


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MOG"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchatMvO7 <- identifyOverExpressedGenes(cellchatMvO7, group.dataset = "datasets", 
                                             pos.dataset = "MOG", features.name = features.name, 
                                             only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
netMvO7 <- netMappingDEG(cellchatMvO7, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in Day10
net.upMvO7 <- subsetCommunication(cellchatMvO7, net = netMvO7, datasets = "MOG",ligand.logFC = 0.2)
net.upMvO7LR <- subsetCommunication(cellchatMvO7, net = netMvO7, datasets = "MOG",ligand.logFC = 0.2, receptor.logFC = 0.2)
#net.down <- subsetCommunication(cellchatScDays, net = net, datasets = "Day6",ligand.logFC = 0.2)


#write.csv(net.up,paste(plotPath,"/net.up_Day7_MogvsOVA.csv", sep = ""), row.names = T)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in Day6, i.e.,downregulated in Day10
net.downMvO7 <- subsetCommunication(cellchatMvO7, net = netMvO7, datasets = "OVA",ligand.logFC = -0.2, receptor.logFC = -0.2)
net.downMvO7L <- subsetCommunication(cellchatMvO7, net = netMvO7, datasets = "OVA",ligand.logFC = -0.2)

#write.csv(net.down,paste(plotPath,"/net.down_Day7_MogvsOVA.csv", sep = ""), row.names = T)

gene.upMvO7 <- extractGeneSubsetFromPair(net.upMvO7, cellchatMvO7)
gene.downMvO7 <- extractGeneSubsetFromPair(net.downMvO7, cellchatMvO7)


pathways.show <- c("CCL") 
png(file = "CCL_Genechorddiagram_Upreg_D7_MogvOVA.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list3[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.upMvO7, #small.gap = 3, big.gap = 1,
                     scale = T,
                     title.name = paste(pathways.show,"pathways upregulated in" ,names(object.list3)[1]),
                     lab.cex = 1)
dev.off()

#png(file = "CCL_Genechorddiagram_Downreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list3[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.downMvO7, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list3)[1]),
                     lab.cex = 1)
dev.off()

pathways.show <- c( "LAMININ")
weight.max <- getMaxWeight(object.list3, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list3)) {
  netVisual_aggregate(object.list3[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list3)[i], unique(cellchatMvO7@meta$Day)))
}


pathways.show <- c( "ITGAL-ITGB2") 
weight.max <- getMaxWeight(object.list3, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list3)) {
  netVisual_aggregate(object.list3[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list3)[i], unique(cellchatMvO7@meta$Day)))
}

dev.off()

pathways.show <- c( "ICAM") 
weight.max <- getMaxWeight(object.list3, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list3)) {
  netVisual_aggregate(object.list3[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list3)[i], unique(cellchatMvO7@meta$Day)))
}


pathways.show <- c( "CCL") 
weight.max <- getMaxWeight(object.list3, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list3)) {
  netVisual_aggregate(object.list3[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list3)[i], unique(cellchatMvO7@meta$Day)))
}


pathways.show <- c( "CCL") 
getwd()
#png(file = "CCL_Genechorddiagram_Up_SC_overtime.png",units = "in",width = 8, height = 7, res = 400)
#pdf(file = "CCL_Genechorddiagram_Up_SC_overtime.pdf",width = 8, height = 7)

netVisual_chord_gene(object.list3[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.upMvO7, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated \nin MOG scaffold at Day 7" ),
                     lab.cex = 1)

dev.off()
getwd()
saveRDS(cellchatMvO7, file = "cellchat_scaf_compDay7MogvOva.rds")


#Load scaffold D7 MOG vs D6 SC ----
#cellchat.D7MOG <- readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/cellchat_D7MOG.rds")
cellchat.D7MOG@meta$Tissue = "Scaffold"
cellchat@meta$Tissue = "Spinal Cord"
levels(cellchat@idents)
levels(cellchat.D7MOG@idents)


# identity = data.frame(cellchat.D7MOG@meta$ident, row.names = rownames(cellchat.D7MOG@meta))
# colnames(identity) = "short.2"
# identity$short.2 = as.factor(identity$short.2)
# identity$short.2
# identity$short.2 <- revalue(identity$short.2, 
#                   c("Comp Mac" = "Mac", 
#                     "NK" = "NKs",
#                     "CD11b+ DC" = "DCs",
#                     "Helper NK" = "NKs",
#                     "B Cell" = "B Cells",
#                     "CD103+ DC" = "DCs",
#                     "Inf Mon" = "Mon",
#                     "Il4i1+ DC" = "DCs",
#                     "Naive CD4" = "CD4 T",
#                     "Treg" = "CD4 T",
#                     "CD4 Tem" = "CD4 T",
#                     "CD8 Tem" = "CD8 T",
#                     "Naive CD8" = "CD8 T"
#                     )
#                   )
# 
# updateClusterLabels(cellchat.D7MOG,
#                     old.cluster.name = levels(cellchat.D7MOG@idents),
#                     new.cluster.name = (c(
#                       "Mon", "Mac", "Mac", "Neut", "NKs", "Stromal", "DCs", "NKs", "B Cells", "DCs",
#                       "Mon", "DCs", "Plasmablast", "CD4 T", "CD4 T", "CD8 T", "CD8 T", "gd", "CD4 T",
#                       "Fer Myeloid"
#                     ))
# )
# 
# length(dimnames(cellchat.D7MOG@net$prob)[[1]])
# length(c(
#   "Mon", "Mac", "Mac", "Neut", "NKs", "Stromal", "DCs", "NKs", "B Cells", "DCs",
#   "Mon", "DCs", "Plasmablast", "CD4 T", "CD4 T", "CD8 T", "CD8 T", "gd", "CD4 T",
#   "Fer Myeloid"
# ))
# 
# levels(cellchat.D7MOG@idents)
# 
# rownames(cellchat@meta)
# 
# 
# identity = data.frame(group = identity$short.2, row.names = names(object$cellgroup))
# unique(identity$group)
# cellchat.D7MOG <- addMeta(cellchat.D7MOG, meta = identity, meta.name = "short.2")
# levels(cellchat.D7MOG@idents) #check idents are correct
# class(cellchat@meta$ident)
# cellchat.D7MOG@meta$short.2
# names(cellchat.D7MOG@meta)
# cellchat.D7MOG <- setIdent(cellchat.D7MOG, ident.use = "short.2")# set "cell_type" as default cell identity
# levels(cellchat.D7MOG@idents)
# cellchat.D7MOG@meta$cell_type
# 


# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. 
# If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents),levels(cellchat.E13@idents))`
group.new = union(levels(cellchat@idents),levels(cellchat.D7MOG@idents))
cellchat1 <- CellChat::liftCellChat(cellchat, group.new)
levels(cellchat1@idents)

cellchat.D7MOG1 <- CellChat::liftCellChat(cellchat.D7MOG, group.new)
levels(cellchat.D7MOG1@idents)
#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...
# Of note, we did not apply `liftCellChat` to cellchat.E14 here because it contains all cell groups. 
object.list2 <- list(`Spinal Cord` = cellchat1, Scaffold = cellchat.D7MOG1)

cellchatEarlyScaf_Sc <- CellChat::mergeCellChat(object.list2, add.names = names(object.list2), cell.prefix = TRUE)


updateClusterLabels(cellchatEarlyScaf_Sc, 
                    new.order = group.new)


#> Warning in mergeCellChat(object.list, add.names = names(object.list),
#> cell.prefix = TRUE): Prefix cell names!
#> The cell barcodes in merged 'meta' is  rep1_AAACCTGCACCAACCG rep1_AAACGGGAGCCGATTT rep1_AAACGGGAGTATCGAA rep1_AAACGGGCATCTCCCA rep1_AAAGATGCACTTGGAT rep1_AAAGATGCAGTTCATG
#> Warning in mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE): The cell barcodes in merged 'meta' is different from those in the used data matrix.
#>               We now simply assign the colnames in the data matrix to the rownames of merged 'mata'!
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

unique(cellchatEarlyScaf_Sc@meta$datasets)
cellchatEarlyScaf_Sc
head(object.list2$`Spinal Cord`@netP)
head(object.list2$Scaffold@netP)


unique(cellchatEarlyScaf_Sc@meta$ident)


gg1 <- compareInteractions(cellchatEarlyScaf_Sc, show.legend = F, group = c(1,2) )
gg2 <- compareInteractions(cellchatEarlyScaf_Sc, show.legend = F, group = c(1,2), measure = "weight",  )
gg1 + gg2

# two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased
# (or decreased) signaling in the second dataset compared to the first one.
#red more in Day 10
# blue less in Day 10

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchatEarlyScaf_Sc, weight.scale = T)
netVisual_diffInteraction(cellchatEarlyScaf_Sc, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchatEarlyScaf_Sc)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchatEarlyScaf_Sc, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list2, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list2)) {
  netVisual_circle(object.list2[[i]]@net$count, weight.scale = T, 
                   label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list2)[i]))
}


num.link <- sapply(object.list2, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list2)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list2[[i]], title = names(object.list2)[i], 
                                               weight.MinMax = weight.MinMax) + 
    ylim(0,10) + xlim(0,15)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- rankNet(cellchatEarlyScaf_Sc, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchatEarlyScaf_Sc, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

pdf(file = "EarlyPathwaysScvScaf.pdf",width = 8, height = 7)
gg1 + gg2
dev.off()

object.list2


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Spinal Cord"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchatEarlyScaf_Sc <- identifyOverExpressedGenes(cellchatEarlyScaf_Sc, group.dataset = "datasets", 
                                           pos.dataset = "Spinal Cord", features.name = features.name, 
                                           only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net_Early_Sc_Scaf <- netMappingDEG(cellchatEarlyScaf_Sc, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in Day10
net.up_Early_Sc_Scaf <- subsetCommunication(cellchatEarlyScaf_Sc, net = net_Early_Sc_Scaf, datasets = "Spinal Cord",
                                            ligand.logFC = 0.2)
#net.down <- subsetCommunication(cellchatEarlyScaf_Sc, net = net, datasets = "Day6",ligand.logFC = 0.2)

unique(net.up_Early_Sc_Scaf$source)

#write.csv(net.up,paste(plotPath,"/net.up_Day7_MogvsOVA.csv", sep = ""), row.names = T)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in Day6, i.e.,downregulated in Day10
net.down_Early_Sc_Scaf <- subsetCommunication(cellchatEarlyScaf_Sc, net = net_Early_Sc_Scaf,
                                              datasets = "Scaffold",ligand.logFC = -0.2, receptor.logFC = -0.2)

#write.csv(net.down,paste(plotPath,"/net.down_Day7_MogvsOVA.csv", sep = ""), row.names = T)


pathways.show <- c("CCL") 
#png(file = "CCL_Genechorddiagram_Upreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list2[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.up_Early_Sc_Scaf, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated in" ,names(object.list2)[1]),
                     lab.cex = 1)
dev.off()

#png(file = "CCL_Genechorddiagram_Downreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list2[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.down_Early_Sc_Scaf, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list2)[1]),
                     lab.cex = 1)
dev.off()

pathways.show <- c( "LAMININ")
weight.max <- getMaxWeight(object.list2, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list2)) {
  netVisual_aggregate(object.list2[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list2)[i]))
}

dev.off()

netVisual_aggregate(object.list2[[2]], signaling = pathways.show, layout = "circle", 
                    edge.weight.max = weight.max[1], edge.width.max = 10, 
                    signaling.name = paste(pathways.show, names(object.list2)[2]))

netVisual_chord_gene(object.list2[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = "LAMININ",legend.pos.x = 8,
                     slot.name = 'net',net = net.upMvO7,
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show ),
                     lab.cex = 1)


pathways.show <- c( "ITGAL-ITGB2") 
weight.max <- getMaxWeight(object.list2, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list2)) {
  netVisual_aggregate(object.list2[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list2)[i], unique(cellchat@meta$Day)))
}

dev.off()

pathways.show <- c( "ICAM") 
weight.max <- getMaxWeight(object.list2, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list2)) {
  netVisual_aggregate(object.list2[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list2)[i]))
}



pathways.show <- c( "CCL") 
group.new
weight.max <- getMaxWeight(object.list2, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list2)) {
  i = 2
  netVisual_aggregate(object.list2[[i]], signaling = pathways.show, layout = "circle", 
                      color.use = c("darkorchid4","#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270",
                                    "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen", "springgreen", 
                                    "bisque", "azure4"),
                      #idents.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                      #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                      #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                      edge.weight.max = weight.max[1], edge.width.max = 5, remove.isolate = F,
                      signaling.name = paste(pathways.show, names(object.list2)[i], unique(cellchat@meta$Day)))
}

netVisual_aggregate(cellchat.D7MOG, signaling = pathways.show, layout = "circle", 
                    color.use = c("darkorchid4","#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270",
                                  "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen", "springgreen", 
                                  "bisque", "azure4"),
                    #idents.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    edge.weight.max = weight.max[1], edge.width.max = 5, remove.isolate = F, 
                    #cell.order = levels(cellchat.D7MOG@idents),
                    signaling.name = paste(pathways.show, names(object.list2)[i], unique(cellchat@meta$Day)))

par(mfrow = c(1,2), xpd=TRUE)
levels(cellchat.D7MOG@idents)
netVisual_aggregate(cellchat.D7MOG, signaling = pathways.show, layout = "circle", 
                    color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270",
                                  "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                    idents.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    edge.weight.max = weight.max[1], edge.width.max = 5, remove.isolate = T, 
                    #cell.order = levels(cellchat.D7MOG@idents),
                    signaling.name = paste(pathways.show, names(object.list2)[i], unique(cellchat@meta$Day)))
dev.off()

par(mfrow = c(1,2), xpd=TRUE)

common_interactionsEarly = intersect(df.net$interaction_ST, cellchat.D7MOG.net$interaction_ST)
intersect(df.net.ccl.cxcl$interaction_ST, Mog_notinOVAD7$interaction_ST)
intersect(common_ccl_MOGD7$interaction_ST, common_ccl_ScD6$interaction_ST)
common_ccl_MOGD7 = cellchat.D7MOG.net.ccl.cxcl[which(cellchat.D7MOG.net.ccl.cxcl$interaction_ST %in% common_interactionsEarly),]
Mog_notinOVAD7 = cellchat.D7MOG.net[-which(cellchat.D7MOG.net$interaction_ST %in% cellchat.D7OVA.net$interaction_ST ),]

length(cellchat.D7OVA.net.ccl.cxcl$interaction_ST)
length(cellchat.D7MOG.net.ccl.cxcl$interaction_ST)
length(cellchat.D7OVA.net.ccl.cxcl$interaction_ST %in% cellchat.D7MOG.net.ccl.cxcl$interaction_ST )
length(cellchat.D7MOG.net.ccl.cxcl$interaction_ST %in% cellchat.D7OVA.net.ccl.cxcl$interaction_ST )

common_ccl_OVAD7 = cellchat.D7OVA.net.ccl.cxcl[which(cellchat.D7OVA.net.ccl.cxcl$interaction_ST %in% common_interactionsEarly),]

intersect(common_ccl_OVAD7$interaction_ST, common_ccl_MOGD7$interaction_ST)
intersect(common_ccl_MOGD7$interaction_ST, common_ccl_OVAD7$interaction_ST)
intersect(Mog_notinOVAD7$interaction_ST, common_ccl_ScD6$interaction_ST)

common_ccl_ScD6 = df.net.ccl.cxcl[which(df.net.ccl.cxcl$interaction_ST %in% common_interactionsEarly),]

common_ccl_ScD6_mog_notOVA = intersect(Mog_notinOVAD7$interaction_ST, common_ccl_ScD6$interaction_ST)
common_ccl_ScD6_mog_notOVA_net = cellchat.D7MOG.net.ccl.cxcl[which(cellchat.D7MOG.net.ccl.cxcl$interaction_ST %in% common_ccl_ScD6_mog_notOVA ),]

ccl_ScD6_notOVA_net = df.net.ccl.cxcl[-which(df.net.ccl.cxcl$interaction_ST %in% cellchat.D7OVA.net.ccl.cxcl$interaction_ST ),]
ccl_ScD6_notOVAoMOG_net = ccl_ScD6_notOVA_net[-which( ccl_ScD6_notOVA_net$interaction_ST %in% common_ccl_MOGD7$interaction_ST),]

intersect(ccl_ScD6_notOVA_net$interaction_ST, cellchat.D7OVA.net.ccl.cxcl$interaction_ST)

#MvOD7.netup = net.up
net.upMvO7$interaction_ST = paste0(net.upMvO7$interaction_name, "_",
                                   net.upMvO7$source, "_",
                                   net.upMvO7$target)

net.downMvO7$interaction_ST = paste0(net.downMvO7$interaction_name, "_",
                                   net.downMvO7$source, "_",
                                   net.downMvO7$target)

common_ccl_ScD6_upmogvOVA = df.net.ccl.cxcl[which( df.net.ccl.cxcl$interaction_ST %in% net.upMvO7$interaction_ST),]

length(common_ccl_ScD6_upmogvOVA$source)
length(net.upMvO7[which(net.upMvO7$interaction_ST %in% df.net.ccl.cxcl$interaction_ST),1])

png(file = "CommonCCL_Genechorddiagram_Scaf_early_upMogvOVA_inSC.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "CommonCCL_Genechorddiagram_Scaf_early_upMogvOVA_inSC.pdf",width = 8, height = 7)

netVisual_chord_gene(object.list2[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = "CCL",legend.pos.x = 8,
                     slot.name = 'net',net = common_ccl_ScD6_upmogvOVA,
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways in\n both the spinal cord and scaffold" ),
                     lab.cex = 1)
dev.off()


Dif_ccl_MOGD7 = cellchat.D7MOG.net.ccl.cxcl[-which(cellchat.D7MOG.net.ccl.cxcl$interaction_ST %in% common_ccl_ScD6_upmogvOVA$interaction_ST),]
Dif_ccl_MOGD7_notOVA = Dif_ccl_MOGD7[-which(Dif_ccl_MOGD7$interaction_ST %in% cellchat.D7OVA.net.ccl.cxcl$interaction_ST),]

Dif_ccl_OVAD7 = cellchat.D7OVA.net.ccl.cxcl[-which(cellchat.D7OVA.net.ccl.cxcl$interaction_ST %in% common_interactionsEarly),]

Dif_ccl_ScD6 = df.net.ccl.cxcl[-which(df.net.ccl.cxcl$interaction_ST %in% cellchat.D7MOG.net.ccl.cxcl$interaction_ST),]
Dif_ccl_ScD61 = Dif_ccl_ScD6[-which(Dif_ccl_ScD6$interaction_ST %in% cellchat.D7OVA.net.ccl.cxcl$interaction_ST),]

pathways.show = "CCL"
png(file = "DifCCL_Genechorddiagram_SC_early_notinOVAoMOG.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "DifCCL_Genechorddiagram_SC_early_notinOVAoMOG.pdf",width = 8, height = 7)

netVisual_chord_gene(object.list2[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = "CCL",legend.pos.x = 8,
                     slot.name = 'net', net = Dif_ccl_ScD61, 
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways only in the spinal cord" ),
                     lab.cex = 1)
dev.off()



png(file = "DifCCL_Genechorddiagram_Scaf_early_notinOVA.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "DifCCL_Genechorddiagram_Scaf_early_notinOVA.pdf",width = 8, height = 7)


netVisual_chord_gene(object.list2[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = c("CCL"),legend.pos.x = 8,
                     slot.name = 'net',net = Dif_ccl_MOGD7_notOVA, 
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways only in the Scaffold" ),
                     lab.cex = 1)
dev.off()



png(file = "CCL_Genechorddiagram_SC_early.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list2[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net', #net = Dif_ccl_ScD6, 
                     small.gap = 3, big.gap = 1, scale = T,
                     sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show,"pathways only in the spinal cord" ),
                     lab.cex = 1)
dev.off()



png(file = "CCL_Genechorddiagram_Scaf_early.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list2[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net', #net = Dif_ccl_MOGD7, 
                     small.gap = 3, big.gap = 1, scale = T,
                     sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show,"pathways only in the Scaffold" ),
                     lab.cex = 1)
dev.off()





#D9 scaf ----
D9 = subset(data, idents = 'D9')


Idents(D9) = D9@meta.data$Antigen
mog9 = subset(D9, idents = "MOG")
ova9 = subset(D9, idents = "OVA")

Idents(mog9) = mog9@meta.data$short
mog9 = SetIdent(mog9, value = mog9@meta.data$short)

#MOG Day 9 ----
cellchat.D9MOG <- createCellChat(object = mog9, group.by = 'short', assay = "RNA")
cell.labels = levels(cellchat.D9MOG@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat.D9MOG@idents)) # number of cells in each cell group
groupSize

# set the used database in the object
cellchat.D9MOG@DB <- CellChatDB.use

cellchat.D9MOG <- subsetData(cellchat.D9MOG) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat.D9MOG <- identifyOverExpressedGenes(cellchat.D9MOG)
cellchat.D9MOG <- identifyOverExpressedInteractions(cellchat.D9MOG)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.mouse, raw.use = F)

# Compute Communicaiton probability and infer network ---------------------
cellchat.D9MOG <- computeCommunProb(cellchat.D9MOG) #remove pop.size = F
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.D9MOG <- filterCommunication(cellchat.D9MOG, min.cells = 4)
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-08-22 22:38:49.873196]"
# |=================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-08-22 22:41:50.0529]"
# Warning message:
#   In UseMethod("depth") :
#   no applicable method for 'depth' applied to an object of class "NULL"


# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
cellchat.D9MOG.net <- subsetCommunication(cellchat.D9MOG) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

cellchat.D9MOG.net$interaction_ST = paste0(cellchat.D9MOG.net$interaction_name, "_",
                                                    cellchat.D9MOG.net$source, "_",
                                                    cellchat.D9MOG.net$target)


cellchat.D9MOG.net.path <- subsetCommunication(cellchat.D9MOG, slot.name = "netP") # look at level of signaling pathways

cellchat.D9MOG.net.ccl.cxcl <- subsetCommunication(cellchat.D9MOG, signaling = c("CCL", "CXCL"))
cellchat.D9MOG.net.ccl.cxcl$interaction_ST = paste0(cellchat.D9MOG.net.ccl.cxcl$interaction_name, "_",
                                                    cellchat.D9MOG.net.ccl.cxcl$source, "_",
                                                    cellchat.D9MOG.net.ccl.cxcl$target)


# Infer cell-cell communication at signaling pathway level ----------------
cellchat.D9MOG <- computeCommunProbPathway(cellchat.D9MOG)


# Calculate aggregated communication network ------------------------------

cellchat.D9MOG <- aggregateNet(cellchat.D9MOG)

groupSize <- as.numeric(table(cellchat.D9MOG@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.D9MOG@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.D9MOG@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat.D9MOG@net$weight
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(3,4), mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


cellchat.D9MOG@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat.D9MOG@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat.D9MOG@idents)

netVisual_chord_gene(cellchat.D9MOG, signaling = c("CCL","CXCL"),legend.pos.x = 2, #sources.use = c(4,5),
                     small.gap = 3, big.gap = 10)
dev.off()

# Part IV: Systems analysis of cell-cell communication network ----
# Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
# (A) Compute and visualize the network centrality scores

# Compute the network centrality scores
cellchat.D9MOG <- netAnalysis_computeCentrality(cellchat.D9MOG, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.D9MOG, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# (B) Visualize dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat.D9MOG)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat.D9MOG, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

#OVA D9 ----
Idents(ova9) = ova9@meta.data$short
ova9 = SetIdent(ova9, value = ova9@meta.data$short)
cellchat.D9OVA <- createCellChat(object = ova9, group.by = 'short', assay = "RNA")
levels(cellchat.D9OVA@idents) # show factor levels of the cell labels
as.numeric(table(cellchat.D9OVA@idents)) # number of cells in each cell group


# set the used database in the object
cellchat.D9OVA@DB <- CellChatDB.use

cellchat.D9OVA <- subsetData(cellchat.D9OVA) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat.D9OVA <- identifyOverExpressedGenes(cellchat.D9OVA)
cellchat.D9OVA <- identifyOverExpressedInteractions(cellchat.D9OVA)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.mouse, raw.use = F)

# Compute Communicaiton probability and infer network ---------------------
cellchat.D9OVA <- computeCommunProb(cellchat.D9OVA) #remove pop.size = F
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.D9OVA <- filterCommunication(cellchat.D9OVA, min.cells = 4)
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-08-22 22:38:49.873196]"
# |=================================================================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-08-22 22:41:50.0529]"
# Warning message:
#   In UseMethod("depth") :
#   no applicable method for 'depth' applied to an object of class "NULL"


# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
cellchat.D9OVA.net <- subsetCommunication(cellchat.D9OVA) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

cellchat.D9OVA.net$interaction_ST = paste0(cellchat.D9OVA.net$interaction_name, "_",
                                                    cellchat.D9OVA.net$source, "_",
                                                    cellchat.D9OVA.net$target)

cellchat.D9OVA.net.path <- subsetCommunication(cellchat.D9OVA, slot.name = "netP") # look at level of signaling pathways

cellchat.D9OVA.net.ccl.cxcl <- subsetCommunication(cellchat.D9OVA, signaling = c("CCL", "CXCL"))
cellchat.D9OVA.net.ccl.cxcl$interaction_ST = paste0(cellchat.D9OVA.net.ccl.cxcl$interaction_name, "_",
                                                    cellchat.D9OVA.net.ccl.cxcl$source, "_",
                                                    cellchat.D9OVA.net.ccl.cxcl$target)


#common_interactionsEarly = intersect(df.net.ccl.cxcl$interaction_ST, cellchat.D7MOG.net.ccl.cxcl$interaction_ST)

# Infer cell-cell communication at signaling pathway level ----------------
cellchat.D9OVA <- computeCommunProbPathway(cellchat.D9OVA)


# Calculate aggregated communication network ------------------------------

cellchat.D9OVA <- aggregateNet(cellchat.D9OVA)

groupSize <- as.numeric(table(cellchat.D9OVA@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.D9OVA@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.D9OVA@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat.D9OVA@net$weight
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(3,4), mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()


cellchat.D9OVA@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat.D9OVA@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat.D9OVA@idents)

netVisual_chord_gene(cellchat.D9OVA, signaling = c("CCL","CXCL"),legend.pos.x = 2, #sources.use = c(4,5),
                     small.gap = 3, big.gap = 10)
dev.off()

# Part IV: Systems analysis of cell-cell communication network ----
# Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
# (A) Compute and visualize the network centrality scores

# Compute the network centrality scores
cellchat.D9OVA <- netAnalysis_computeCentrality(cellchat.D9OVA, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.D9OVA, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# (B) Visualize dominant senders (sources) and receivers (targets) in a 2D space

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat.D9OVA)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat.D9OVA, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

saveRDS(cellchat.D9OVA, file = "cellchat_D9OVA.rds")
saveRDS(cellchat.D9MOG, file = "cellchat_D9MOG.rds")
saveRDS(cellchat.D7OVA, file = "cellchat_D7OVA.rds")
saveRDS(cellchat.D7MOG, file = "cellchat_D7MOG.rds")

#Compare Day 9 MOG v OVA ----

object.list4 <- list(MOG = cellchat.D9MOG, OVA = cellchat.D9OVA)
cellchatMvO9 <- mergeCellChat(object.list4, add.names = names(object.list4))

unique(cellchatMvO9@meta)

gg1 <- compareInteractions(cellchatMvO9, show.legend = F, group = c(1,2) )
gg2 <- compareInteractions(cellchatMvO9, show.legend = F, group = c(1,2), measure = "weight" )
gg1 + gg2

# two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased
# (or decreased) signaling in the second dataset compared to the first one.
#red more in Day 10
# blue less in Day 10

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchatMvO9, weight.scale = T)
netVisual_diffInteraction(cellchatMvO9, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchatMvO9)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchatMvO9, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list4, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list4)) {
  netVisual_circle(object.list4[[i]]@net$count, weight.scale = T, 
                   label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list4)[i]))
}


num.link <- sapply(object.list4, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list4)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list4[[i]], 
                                               title = names(object.list4)[i],
                                               weight.MinMax = weight.MinMax) + xlim(0,12) + ylim(0,10)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- rankNet(cellchatMvO9, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchatMvO9, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

object.list4


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MOG"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchatMvO9 <- identifyOverExpressedGenes(cellchatMvO9, group.dataset = "datasets", 
                                           pos.dataset = "MOG", features.name = features.name, 
                                           only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
netMvO9 <- netMappingDEG(cellchatMvO9, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in Day10
net.upMvO9 <- subsetCommunication(cellchatMvO9, net = netMvO9, datasets = "MOG",ligand.logFC = 0.2)
#net.down <- subsetCommunication(cellchatScDays, net = net, datasets = "Day6",ligand.logFC = 0.2)


#write.csv(net.up,paste(plotPath,"/net.up_Day7_MogvsOVA.csv", sep = ""), row.names = T)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in Day6, i.e.,downregulated in Day10
net.downMvO9 <- subsetCommunication(cellchatMvO9, net = netMvO9, datasets = "OVA",ligand.logFC = -0.2, receptor.logFC = -0.2)

#write.csv(net.down,paste(plotPath,"/net.down_Day7_MogvsOVA.csv", sep = ""), row.names = T)


pathways.show <- c("CCL") 
png(file = "CCL_Genechorddiagram_Upreg_D9_MogvOVA.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list4[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.upMvO9, #small.gap = 3, big.gap = 1,
                     scale = T,
                     title.name = paste(pathways.show,"pathways upregulated in" ,names(object.list4)[1]),
                     lab.cex = 1)
dev.off()

#png(file = "CCL_Genechorddiagram_Downreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list4[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.downMvO9, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list4)[1]),
                     lab.cex = 1)
dev.off()

pathways.show <- c( "LAMININ")
weight.max <- getMaxWeight(object.list4, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list4)) {
  netVisual_aggregate(object.list4[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list4)[i], unique(cellchatMvO9@meta$Day)))
}


pathways.show <- c( "ITGAL-ITGB2") 
weight.max <- getMaxWeight(object.list4, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list4)) {
  netVisual_aggregate(object.list4[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list4)[i], unique(cellchatMvO9@meta$Day)))
}

dev.off()

pathways.show <- c( "ICAM") 
weight.max <- getMaxWeight(object.list4, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list4)) {
  netVisual_aggregate(object.list4[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list4)[i], unique(cellchatMvO9@meta$Day)))
}


pathways.show <- c( "CCL") 
weight.max <- getMaxWeight(object.list4, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list4)) {
  netVisual_aggregate(object.list4[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list4)[i], unique(cellchatMvO9@meta$Day)))
}

pathways.show <- c( "CCL") 
getwd()
#png(file = "CCL_Genechorddiagram_Up_SC_overtime.png",units = "in",width = 8, height = 7, res = 400)
#pdf(file = "CCL_Genechorddiagram_Up_SC_overtime.pdf",width = 8, height = 7)

netVisual_chord_gene(object.list4[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.upMvO9, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated in MOG" ),
                     lab.cex = 1)

dev.off()
saveRDS(cellchatMvO9, file = "cellchat_scaf_compDay9MogvOva.rds")



#Compare Day 7 v 9  MOG ----

object.list6 <- list(D7 = cellchat.D7MOG, D9 = cellchat.D9MOG)
cellchatM7v9 <- mergeCellChat(object.list6, add.names = names(object.list6))

unique(cellchatM7v9@meta)

gg1 <- compareInteractions(cellchatM7v9, show.legend = F, group = c(1,2) )
gg2 <- compareInteractions(cellchatM7v9, show.legend = F, group = c(1,2), measure = "weight" )
gg1 + gg2

# two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased
# (or decreased) signaling in the second dataset compared to the first one.
#red more in Day 10
# blue less in Day 10

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchatM7v9, weight.scale = T)
netVisual_diffInteraction(cellchatM7v9, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchatM7v9)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchatM7v9, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list6, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list6)) {
  netVisual_circle(object.list6[[i]]@net$count, weight.scale = T, 
                   label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list6)[i]))
}


num.link <- sapply(object.list6, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list6)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list6[[i]], 
                                               title = names(object.list6)[i],
                                               weight.MinMax = weight.MinMax) + xlim(0,12) + ylim(0,10)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- rankNet(cellchatM7v9, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchatM7v9, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

object.list6


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "D9"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchatM7v9 <- identifyOverExpressedGenes(cellchatM7v9, group.dataset = "datasets", 
                                           pos.dataset = "D9", features.name = features.name, 
                                           only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
netM7v9 <- netMappingDEG(cellchatM7v9, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in Day10
net.upM7v9 <- subsetCommunication(cellchatM7v9, net = netM7v9, datasets = "D9",ligand.logFC = 0.2)
#net.down <- subsetCommunication(cellchatScDays, net = net, datasets = "Day6",ligand.logFC = 0.2)


#write.csv(net.up,paste(plotPath,"/net.up_Day7_MogvsOVA.csv", sep = ""), row.names = T)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in Day6, i.e.,downregulated in Day10
net.downM7v9 <- subsetCommunication(cellchatM7v9, net = netM7v9, datasets = "D7",ligand.logFC = -0.2, receptor.logFC = -0.2)

#write.csv(net.down,paste(plotPath,"/net.down_Day7_MogvsOVA.csv", sep = ""), row.names = T)


pathways.show <- c("CCL") 
#png(file = "CCL_Genechorddiagram_Upreg_D9_MogvOVA.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list6[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.upM7v9, #small.gap = 3, big.gap = 1,
                     scale = T,
                     title.name = paste(pathways.show,"pathways upregulated in" ,names(object.list6)[2]),
                     lab.cex = 1)
dev.off()

#png(file = "CCL_Genechorddiagram_Downreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list6[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.downM7v9, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list6)[2]),
                     lab.cex = 1)


netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.down, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list)[2]),
                     lab.cex = 1)
dev.off()

pathways.show <- c( "LAMININ")
weight.max <- getMaxWeight(object.list6, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list6)) {
  netVisual_aggregate(object.list6[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list6)[i]))
}


pathways.show <- c( "ITGAL-ITGB2") 
weight.max <- getMaxWeight(object.list6, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list6)) {
  netVisual_aggregate(object.list6[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list6)[i]))
}


dev.off()

pathways.show <- c( "ICAM") 
weight.max <- getMaxWeight(object.list6, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list6)) {
  netVisual_aggregate(object.list6[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list6)[i]))
}


pathways.show <- c( "CCL") 
weight.max <- getMaxWeight(object.list6, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list6)) {
  netVisual_aggregate(object.list6[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list6)[i]))
}

pathways.show <- c( "CCL") 
getwd()
#png(file = "CCL_Genechorddiagram_Up_SC_overtime.png",units = "in",width = 8, height = 7, res = 400)
#pdf(file = "CCL_Genechorddiagram_Up_SC_overtime.pdf",width = 8, height = 7)

netVisual_chord_gene(object.list6[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.upM7v9, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated at D9" ),
                     lab.cex = 1)

dev.off()
saveRDS(cellchatM7v9, file = "cellchat_scaf_compDay9v7Mog.rds")

#Load scaffold D9 MOG vs 10 SC ----
#cellchat.D7MOG <- readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/cellchat_D7MOG.rds")
cellchat.D9MOG@meta$Tissue = "Scaffold"
cellchat10@meta$Tissue = "Spinal Cord"
levels(cellchat10@idents)
levels(cellchat.D9MOG@idents)


# Define the cell labels to lift up by combining both cell labels from the conditions
# Please note that the order of cell groups in `group.new` will affect the appearance order when visualizing the cell-cell communication. 
# If there are unique cell groups in both the conditions, you should define `group.new = union(levels(cellchat.E14@idents),levels(cellchat.E13@idents))`
group.new1 = union(levels(cellchat10@idents),levels(cellchat.D9MOG@idents))
cellchat10_1 <- CellChat::liftCellChat(cellchat10, group.new1)
levels(cellchat10_1@idents)
levels(cellchat10_1@meta$ident)
levels(cellchat1@meta$ident)

cellchat.D9MOG1 <- CellChat::liftCellChat(cellchat.D9MOG, group.new1)
levels(cellchat.D9MOG1@idents)
levels(cellchat.D7MOG1@idents)
levels(cellchat.D9MOG1@meta$ident)


#> The CellChat object will be lifted up using the cell labels FIB-A, FIB-B, FIB-P, DC, Pericyte, MYL, Immune, ENDO, Muscle, MELA, Basal-P, Basal, Spinious
#> Update slots object@net, object@netP, object@idents in a single dataset...
# Of note, we did not apply `liftCellChat` to cellchat.E14 here because it contains all cell groups. 
object.list5 <- list(`Spinal Cord` = cellchat10_1, Scaffold = cellchat.D9MOG1)


levels(object.list5$Scaffold@idents)
levels(cellchat.D9MOG@idents)

unique(object.list5$`Spinal Cord`@idents)
unique(object.list5$`Spinal Cord`@meta$ident)
levels(cellchat10@idents)

object.list5$Scaffold@meta$ident = object.list5$Scaffold@idents
object.list5$`Spinal Cord`@meta$ident = object.list5$`Spinal Cord`@idents


object.list5$`Spinal Cord`@meta$ident[1:10]
object.list5$`Spinal Cord`@idents[1:10]

object.list5$Scaffold@idents[1:10]

object.list5$Scaffold@meta$ident[1:10]

cellchatLateScaf_Sc <- CellChat::mergeCellChat(object.list5, add.names = names(object.list5), cell.prefix = TRUE)

cellchatEarlyScaf_Sc@idents
cellchatLateScaf_Sc@idents

#> Warning in mergeCellChat(object.list, add.names = names(object.list),
#> cell.prefix = TRUE): Prefix cell names!
#> The cell barcodes in merged 'meta' is  rep1_AAACCTGCACCAACCG rep1_AAACGGGAGCCGATTT rep1_AAACGGGAGTATCGAA rep1_AAACGGGCATCTCCCA rep1_AAAGATGCACTTGGAT rep1_AAAGATGCAGTTCATG
#> Warning in mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE): The cell barcodes in merged 'meta' is different from those in the used data matrix.
#>               We now simply assign the colnames in the data matrix to the rownames of merged 'mata'!
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

unique(cellchatLateScaf_Sc@meta$datasets)
cellchatLateScaf_Sc
head(object.list5$`Spinal Cord`@netP)
head(object.list5$Scaffold@netP)


unique(cellchatLateScaf_Sc@meta$ident)


gg1 <- compareInteractions(cellchatLateScaf_Sc, show.legend = F, group = c(1,2) )
gg2 <- compareInteractions(cellchatLateScaf_Sc, show.legend = F, group = c(1,2), measure = "weight",  )
gg1 + gg2

# two datasets can be visualized using circle plot, where red (or blue) colored edges represent increased
# (or decreased) signaling in the second dataset compared to the first one.
#red more in Day 10
# blue less in Day 10

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchatLateScaf_Sc, weight.scale = T)
netVisual_diffInteraction(cellchatLateScaf_Sc, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchatLateScaf_Sc)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchatLateScaf_Sc, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list5, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list5)) {
  netVisual_circle(object.list5[[i]]@net$count, weight.scale = T, 
                   label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(object.list5)[i]))
}


num.link <- sapply(object.list5, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list5)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list5[[i]], title = names(object.list5)[i], 
                                               weight.MinMax = weight.MinMax) + 
    ylim(0,10) + xlim(0,15)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)



num.link <- sapply(object.list5, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list5)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list5[[i]], title = names(object.list5)[i], 
                                               weight.MinMax = weight.MinMax, signaling = "CCL") + 
    ylim(0,2) + xlim(0,2)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- rankNet(cellchatLateScaf_Sc, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchatLateScaf_Sc, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)

gg1 + gg2

object.list5


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Spinal Cord"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchatLateScaf_Sc <- identifyOverExpressedGenes(cellchatLateScaf_Sc, group.dataset = "datasets", 
                                                   pos.dataset = "Spinal Cord", features.name = features.name, 
                                                   only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
netLateScaf_Sc <- netMappingDEG(cellchatLateScaf_Sc, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in Day10
net.upLateScaf_Sc <- subsetCommunication(cellchatLateScaf_Sc, net = netLateScaf_Sc, datasets = "Spinal Cord",ligand.logFC = 0.2)
#net.down <- subsetCommunication(cellchatEarlyScaf_Sc, net = net, datasets = "Day6",ligand.logFC = 0.2)


#write.csv(net.up,paste(plotPath,"/net.up_Day7_MogvsOVA.csv", sep = ""), row.names = T)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in Day6, i.e.,downregulated in Day10
net.downLateScaf_Sc <- subsetCommunication(cellchatLateScaf_Sc, net = netLateScaf_Sc, datasets = "Scaffold",ligand.logFC = -0.2, receptor.logFC = -0.2)

#write.csv(net.down,paste(plotPath,"/net.down_Day7_MogvsOVA.csv", sep = ""), row.names = T)


pathways.show <- c("CCL") 

#png(file = "CCL_Genechorddiagram_Downreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list5[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.downLateScaf_Sc, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list5)[1]),
                     lab.cex = 1)


netVisual_chord_gene(object.list5[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.upLateScaf_Sc, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list5)[1]),
                     lab.cex = 1)
dev.off()

pathways.show <- c( "LAMININ")
weight.max <- getMaxWeight(object.list5, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list5)) {
  netVisual_aggregate(object.list5[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list5)[i], "Late"))
}

netVisual_aggregate(cellchat.D7MOG, signaling = pathways.show, layout = "circle", 
                    edge.weight.max = weight.max[1], edge.width.max = 10, 
                    signaling.name = paste(pathways.show, names(object.list5)[i], "Late"))

netVisual_aggregate(cellchat.D9MOG, signaling = pathways.show, layout = "circle", 
                    edge.weight.max = weight.max[1], edge.width.max = 10, 
                    signaling.name = paste(pathways.show, names(object.list5)[i], "Late"))

netVisual_aggregate(cellchat10, signaling = pathways.show, layout = "circle", 
                    edge.weight.max = weight.max[1], edge.width.max = 10, 
                    signaling.name = paste(pathways.show, unique(cellchat10@meta$day_group)))

netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", 
                    edge.weight.max = weight.max[1], edge.width.max = 10, 
                    signaling.name = paste(pathways.show, unique(cellchat@meta$day_group)))


par(mfrow = c(1,2), xpd=TRUE)
i = 2
netVisual_aggregate(object.list5[[i]], signaling = pathways.show, layout = "circle", 
                    color.use = c("darkorchid4","#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270",
                                  "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen", "springgreen", 
                                  "bisque", "azure4"),
                    #idents.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    edge.weight.max = weight.max[1], edge.width.max = 5, remove.isolate = F,
                    signaling.name = paste(pathways.show, names(object.list5)[i], unique(cellchat@meta$Day)))



netVisual_aggregate(cellchat.D9MOG, signaling = pathways.show, layout = "circle", 
                    color.use = c("darkorchid4","#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270",
                                  "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen", "springgreen", 
                                  "bisque", "azure4"),
                    #idents.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    edge.weight.max = weight.max[1], edge.width.max = 5, remove.isolate = F, 
                    #cell.order = levels(cellchat.D7MOG@idents),
                    signaling.name = paste(pathways.show, names(object.list2)[i], unique(cellchat@meta$Day)))

plotGeneExpression(cellchat.D9MOG1, signaling = "LAMININ", enriched.only = TRUE, type = "violin")


levels(cellchat.D9MOG1@idents)
levels(cellchat.D9MOG@idents)

pathways.show <- c( "ITGAL-ITGB2") 
weight.max <- getMaxWeight(object.list5, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list5)) {
  netVisual_aggregate(object.list5[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list5)[i]))
}

dev.off()

pathways.show <- c( "ICAM") 
weight.max <- getMaxWeight(object.list5, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list5)) {
  netVisual_aggregate(object.list5[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list5)[i]))
}

object.list2[[1]]@netP

pathways.show <- c( "CCL") 
group.new1
weight.max <- getMaxWeight(object.list5, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list5)) {
  netVisual_aggregate(object.list5[[i]], signaling = pathways.show, layout = "circle", 
                      color.use = c("darkorchid4","#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270",
                                    "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen", "springgreen", 
                                    "bisque", "azure4"),
                      #idents.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                      #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                      #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                      edge.weight.max = weight.max[1], edge.width.max = 5, remove.isolate = F,
                      signaling.name = paste(pathways.show, names(object.list5)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
levels(cellchat.D9MOG@idents)
netVisual_aggregate(cellchat.D9MOG, signaling = pathways.show, layout = "circle", 
                    color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270",
                                  "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                    idents.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                    edge.weight.max = weight.max[1], edge.width.max = 5, remove.isolate = T, 
                    #cell.order = levels(cellchat.D7MOG@idents),
                    signaling.name = paste(pathways.show, names(object.list2)[i], unique(cellchat@meta$Day)))
dev.off()

par(mfrow = c(1,2), xpd=TRUE)

#MvOD9.netup = net.up
net.upMvO9$interaction_ST = paste0(net.upMvO9$interaction_name, "_",
                                   net.upMvO9$source, "_",
                                   net.upMvO9$target)

df10.net.ccl.cxcl$interaction_ST = paste0(df10.net.ccl.cxcl$interaction_name, "_",
                                          df10.net.ccl.cxcl$source, "_",
                                          df10.net.ccl.cxcl$target)

common_ccl_ScD10_upmogvOVA = df10.net.ccl.cxcl[which( df10.net.ccl.cxcl$interaction_ST %in% net.upMvO9$interaction_ST),]
length(common_ccl_ScD10_upmogvOVA$source)
length(net.upMvO9[which(net.upMvO9$interaction_ST %in% df10.net.ccl.cxcl$interaction_ST),1])

png(file = "CommonCCL_Genechorddiagram_Scaf_late_upMogvOVA_inSC.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "CommonCCL_Genechorddiagram_Scaf_late_upMogvOVA_inSC.pdf",width = 8, height = 7)
netVisual_chord_gene(object.list5[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = common_ccl_ScD10_upmogvOVA,
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show,"pathways in\n both the spinal cord and scaffold" ),
                     lab.cex = 1)
dev.off()


Dif_ccl_MOGD9 = cellchat.D9MOG.net.ccl.cxcl[-which(cellchat.D9MOG.net.ccl.cxcl$interaction_ST %in% common_ccl_ScD10_upmogvOVA$interaction_ST),]
Dif_ccl_MOGD9_notOVA = cellchat.D9MOG.net.ccl.cxcl[-which(cellchat.D9MOG.net.ccl.cxcl$interaction_ST %in% cellchat.D9OVA.net.ccl.cxcl$interaction_ST),]

#Dif_ccl_OVAD7 = cellchat.D7OVA.net.ccl.cxcl[-which(cellchat.D7OVA.net.ccl.cxcl$interaction_ST %in% common_interactionsEarly),]

Dif_ccl_ScD10 = df10.net.ccl.cxcl[-which(df10.net.ccl.cxcl$interaction_ST %in% cellchat.D9MOG.net.ccl.cxcl$interaction_ST),]
Dif_ccl_ScD10_1 = df10.net.ccl.cxcl[-which(df10.net.ccl.cxcl$interaction_ST %in% cellchat.D9OVA.net.ccl.cxcl$interaction_ST),]


png(file = "DifCCL_Genechorddiagram_SC_late_notinOVAoMOG.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "DifCCL_Genechorddiagram_SC_late_notinOVAoMOG.pdf",width = 8, height = 7)
netVisual_chord_gene(object.list5[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net', net = Dif_ccl_ScD10_1, 
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show,"pathways only in the spinal cord" ),
                     lab.cex = 1)
dev.off()



png(file = "DifCCL_Genechorddiagram_Scaf_late_notinOVA.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "DifCCL_Genechorddiagram_Scaf_late_notinOVA.pdf",width = 8, height = 7)
netVisual_chord_gene(object.list5[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = c("CCL"),legend.pos.x = 8,
                     slot.name = 'net',net = Dif_ccl_MOGD9_notOVA, 
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show,"pathways only in the Scaffold" ),
                     lab.cex = 1)
dev.off()



#png(file = "CCL_Genechorddiagram_SC_early.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list5[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net', #net = Dif_ccl_ScD6, 
                     small.gap = 3, big.gap = 1, scale = T,
                     sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show,"pathways only in the spinal cord" ),
                     lab.cex = 1)
dev.off()



#png(file = "CCL_Genechorddiagram_Scaf_early.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list5[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net', #net = Dif_ccl_MOGD7, 
                     small.gap = 3, big.gap = 1, scale = T,
                     sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show,"pathways only in the Scaffold" ),
                     lab.cex = 1)
dev.off()



#png(file = "CCL_Genechorddiagram_SC_early.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_cell(object.list5[[1]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'netP', #net = net.up, 
                     small.gap = 3, big.gap = 1, scale = T, 
                     sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show,"pathways \nin the spinal cord" ),
                     lab.cex = 1)
dev.off()
#png(file = "CCL_Genechorddiagram_Scaf_early.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_cell(object.list5[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'netP',#net = net.down, 
                     small.gap = 3, big.gap = 1, scale  = T,
                     sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste(pathways.show,"pathways \nin the Scaffold" ),
                     lab.cex = 1)
dev.off()

unique(common.net.up$pathway_name)
pathways.show <- c( "CCL") 
getwd()
#png(file = "CCL_Genechorddiagram_Up_SC_overtime.png",units = "in",width = 8, height = 7, res = 400)
#pdf(file = "CCL_Genechorddiagram_Up_SC_overtime.pdf",width = 8, height = 7)

netVisual_chord_gene(object.list2[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270",
                                   "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated \nin the spinal cord over time" ),
                     lab.cex = 1)
dev.off()


#comparison by individual cellchat objects----
cellchat@meta$Day = cellchat@meta$day_group
cellchat10@meta$Day = cellchat10@meta$day_group
all_data_list = c(cellchat.D7MOG, cellchat.D9MOG, cellchat, cellchat10)

pathways.show <- c( "LAMININ")


weight.max <- getMaxWeight(all_data_list, 
                           slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets


png(file = "Laminin_dotplot_SC_Scaf.png",units = "in",width = 14, height = 14, res = 400)

par(mfrow = c(2,2), xpd=TRUE)

for (i in all_data_list) {
  netVisual_aggregate(i, signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, unique(i@meta$Tissue) ,unique(i@meta$Day)))
}

dev.off()


pathways.show <- c( "ICAM")


weight.max <- getMaxWeight(all_data_list, 
                           slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets


png(file = "ICAM_dotplot_SC_Scaf.png",units = "in",width = 14, height = 14, res = 400)

par(mfrow = c(2,2), xpd=TRUE)

for (i in all_data_list) {
  netVisual_aggregate(i, signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, unique(i@meta$Tissue) ,unique(i@meta$Day)))
}

dev.off()




pathways.show <- c( "ITGAL-ITGB2")


weight.max <- getMaxWeight(all_data_list, 
                           slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets


png(file = "ITGAL_ITGB2_dotplot_SC_Scaf.png",units = "in",width = 14, height = 14, res = 400)

par(mfrow = c(2,2), xpd=TRUE)

for (i in all_data_list) {
  netVisual_aggregate(i, signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, unique(i@meta$Tissue) ,unique(i@meta$Day)))
}

dev.off()


pathways.show <- c( "CCL")


weight.max <- getMaxWeight(all_data_list, 
                           slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets


png(file = "CCL_dotplot_SC_Scaf.png",units = "in",width = 14, height = 14, res = 400)

par(mfrow = c(2,2), xpd=TRUE)

for (i in all_data_list) {
  netVisual_aggregate(i, signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, unique(i@meta$Tissue) ,unique(i@meta$Day)))
}

dev.off()

#comparison of scaf at early and SC at early ----

#df.net = significant interactions at Day 6 SC
#cellchat.D7MOG.net = significant interactions at Day 7 Scaf MOG

#all interactions that are the same Day 6 SC and Day 7 MOG Scaf
common_interactionsEarly = intersect(df.net$interaction_ST, cellchat.D7MOG.net$interaction_ST)
length(intersect(df.net$interaction_ST,net.upMvO7$interaction_ST))
length(intersect(df10.net$interaction_ST,net.upMvO9$interaction_ST))

length(common_interactionsEarly)
length(common_interactionsLate)
common_interactionsEarly.df = df.net[which(df.net$interaction_ST %in% common_interactionsEarly),]
view(common_interactionsEarly.df)
intersect(df.net.ccl.cxcl$interaction_ST, Mog_notinOVAD7$interaction_ST)
intersect(common_ccl_MOGD7$interaction_ST, common_ccl_ScD6$interaction_ST)

#all CCL interactions that are the same Day 6 SC and Day 7 MOG Scaf
common_ccl_MOGD7 = cellchat.D7MOG.net.ccl.cxcl[which(cellchat.D7MOG.net.ccl.cxcl$interaction_ST %in% common_interactionsEarly),]

#all interactions in Day 7 MOG Scaf that are absent in Day 7 OVA Scaf
Mog_notinOVAD7 = cellchat.D7MOG.net[-which(cellchat.D7MOG.net$interaction_ST %in% cellchat.D7OVA.net$interaction_ST ),]

length(cellchat.D7OVA.net.ccl.cxcl$interaction_ST)
length(cellchat.D7MOG.net.ccl.cxcl$interaction_ST)
length(cellchat.D7OVA.net.ccl.cxcl$interaction_ST %in% cellchat.D7MOG.net.ccl.cxcl$interaction_ST )
length(cellchat.D7MOG.net.ccl.cxcl$interaction_ST %in% cellchat.D7OVA.net.ccl.cxcl$interaction_ST )

#all interactions in Day 7 OVA Scaf that are absent in Day 7 OVA Scaf
common_ccl_OVAD7 = cellchat.D7OVA.net.ccl.cxcl[which(cellchat.D7OVA.net.ccl.cxcl$interaction_ST %in% common_interactionsEarly),]

intersect(common_ccl_OVAD7$interaction_ST, common_ccl_MOGD7$interaction_ST)
intersect(common_ccl_MOGD7$interaction_ST, common_ccl_OVAD7$interaction_ST)
intersect(Mog_notinOVAD7$interaction_ST, common_ccl_ScD6$interaction_ST)

#all CCL interactions that are the same Day 6 SC and Day 7 MOG Scaf
common_ccl_ScD6 = df.net.ccl.cxcl[which(df.net.ccl.cxcl$interaction_ST %in% common_interactionsEarly),]

intersect(common_ccl_MOGD7$interaction_ST, common_ccl_ScD6$interaction_ST)
length(common_ccl_MOGD7$interaction_ST)
length(common_ccl_ScD6$interaction_ST)

# CCL intereaactions in MOG and NOT OVA Day 7 scaf, and in Day 6 SC
common_ccl_ScD6_mog_notOVA = intersect(Mog_notinOVAD7$interaction_ST, common_ccl_ScD6$interaction_ST)
common_ccl_ScD6_mog_notOVA_net = cellchat.D7MOG.net.ccl.cxcl[which(cellchat.D7MOG.net.ccl.cxcl$interaction_ST %in% common_ccl_ScD6_mog_notOVA ),]




intersect(ccl_ScD6_notOVA_net$interaction_ST, cellchat.D7OVA.net.ccl.cxcl$interaction_ST)

#MvOD7.netup = net.up
net.upMvO7$interaction_ST = paste0(net.upMvO7$interaction_name, "_",
                                   net.upMvO7$source, "_",
                                   net.upMvO7$target)

net.downMvO7$interaction_ST = paste0(net.downMvO7$interaction_name, "_",
                                     net.downMvO7$source, "_",
                                     net.downMvO7$target)

#CCL interactions in Day 6 SC and upregualted in diseased MOG day 7 Scaf (compared to Day 7 OVA)
common_ccl_ScD6_upmogvOVA = df.net.ccl.cxcl[which( df.net.ccl.cxcl$interaction_ST %in% net.upMvO7$interaction_ST),]

intersect(common_ccl_ScD6_upmogvOVA$interaction_ST, common_ccl_ScD6_mog_notOVA_net$interaction_ST)

#Figure shows CCL interactions in Day 6 SC and upregualted in diseased MOG day 7 Scaf (compared to Day 7 OVA)
png(file = "CommonCCL_Genechorddiagram_Scaf_early_upMogvOVA_inSC_colorcorrect.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "CommonCCL_Genechorddiagram_Scaf_early_upMogvOVA_inSC_colorcorrect.pdf",width = 8, height = 7)

netVisual_chord_gene(cellchat,# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = "CCL",legend.pos.x = 8,
                     slot.name = 'net',net = common_ccl_ScD6_upmogvOVA,
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways in\n both the spinal cord and scaffold" ),
                     lab.cex = 1)
dev.off()

# CCL Pathways in MOG day 7 not in Day 6 SC
Dif_ccl_MOGD7 = cellchat.D7MOG.net.ccl.cxcl[-which(cellchat.D7MOG.net.ccl.cxcl$interaction_ST %in% common_ccl_ScD6_upmogvOVA$interaction_ST),]

# CCL Pathways in MOG day 7 not in Day 6 SC and absent in Day 7 OVA Scaf
Dif_ccl_MOGD7_notOVA = Dif_ccl_MOGD7[-which(Dif_ccl_MOGD7$interaction_ST %in% cellchat.D7OVA.net.ccl.cxcl$interaction_ST),]

# CCL Pathways in MOG day 7 not in Day 6 SC and not upregulated in Day 7 OVA Scaf
Dif_ccl_MOGD7_notupOVA = Dif_ccl_MOGD7[-which(Dif_ccl_MOGD7$interaction_ST %in% net.downMvO7$interaction_ST),]

# CCL Pathways in OVA and not in MOG scaf day 7 or  SC day 6
Dif_ccl_OVAD7 = cellchat.D7OVA.net.ccl.cxcl[-which(cellchat.D7OVA.net.ccl.cxcl$interaction_ST %in% common_interactionsEarly),]

#CCL pathways in SC day 6 not in MOG scaf day 7
Dif_ccl_ScD6 = df.net.ccl.cxcl[-which(df.net.ccl.cxcl$interaction_ST %in% cellchat.D7MOG.net.ccl.cxcl$interaction_ST),]
#CCL pathways in SC day 6 not in MOG scaf day 7 or OVA scaf day 7
Dif_ccl_ScD61 = Dif_ccl_ScD6[-which(Dif_ccl_ScD6$interaction_ST %in% cellchat.D7OVA.net.ccl.cxcl$interaction_ST),]

#CCL pathways in SC day 6 not up in MOG scaf day 7 or down in OVA scaf day 7
Dif_ccl_ScD62 = df.net.ccl.cxcl[-which(df.net.ccl.cxcl$interaction_ST %in% net.upMvO7$interaction_ST),]
#Dif_ccl_ScD62 = Dif_ccl_ScD62[-which(Dif_ccl_ScD62$interaction_ST %in% net.downMvO7$interaction_ST),]

length(Dif_ccl_ScD62$interaction_ST)
length(Dif_ccl_ScD61$interaction_ST)

#figure shows CCL pathways in SC day 6 not in MOG scaf day 7 or OVA scaf day 7
pathways.show = "CCL"
png(file = "DifCCL_Genechorddiagram_SC_early_notinOVAoMOG_colorcorrect.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "DifCCL_Genechorddiagram_SC_early_notinOVAoMOG_colorcorrect.pdf",width = 8, height = 7)

netVisual_chord_gene(cellchat,# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = "CCL",legend.pos.x = 8,
                     slot.name = 'net', net = Dif_ccl_ScD61, 
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways only in the spinal cord" ),
                     lab.cex = 1)
dev.off()


# CCL Pathways in MOG day 7 not in Day 6 SC and not upregulated in Day 7 OVA Scaf
png(file = "DifCCL_Genechorddiagram_Scaf_early_notupinOVA.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "DifCCL_Genechorddiagram_Scaf_early_notupinOVA.pdf",width = 8, height = 7)


netVisual_chord_gene(cellchat.D7MOG,# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = c("CCL"),legend.pos.x = 8,
                     slot.name = 'net',net = Dif_ccl_MOGD7_notupOVA, 
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways only in the Scaffold" ),
                     lab.cex = 1)
dev.off()


# CCL Pathways in MOG day 7 not in Day 6 SC and not upregulated in Day 7 OVA Scaf
png(file = "DifCCL_Genechorddiagram_Scaf_early_notinOVA_colorcorrect.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "DifCCL_Genechorddiagram_Scaf_early_notinOVA_colorcorrect.pdf",width = 8, height = 7)


netVisual_chord_gene(cellchat.D7MOG,# sources.use = c(6) , targets.use = c(Tcells), 
                     color.use = c("#377eb8", "#4daf4a", "#e4211c", "#bc9dcc", "grey", "#f781bf", "blue","#9C5B32", "#E4993F", "#54b0e4",
                                   "#b2df8a","#b2df8a","#e3be00","#e3be00","darkgreen"),
                     signaling = c("CCL"),legend.pos.x = 8,
                     slot.name = 'net',net = Dif_ccl_MOGD7_notOVA, 
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways only in the Scaffold" ),
                     lab.cex = 1)
dev.off()






#comparison of scaf and SC at late ----


#df10.net = significant interactions at Day 10 SC
#cellchat.D9MOG.net = significant interactions at Day 9 Scaf MOG

#all interactions that are the same Day 10 SC and Day 9 MOG Scaf
common_interactionsLate = intersect(df10.net$interaction_ST, cellchat.D9MOG.net$interaction_ST)


#all CCL interactions that are the same Day 10 SC and Day 9 MOG Scaf
common_ccl_MOGD9 = cellchat.D9MOG.net.ccl.cxcl[which(cellchat.D9MOG.net.ccl.cxcl$interaction_ST %in% common_interactionsLate),]

#all interactions in Day 9 MOG Scaf that are absent in Day 9 OVA Scaf
Mog_notinOVAD9 = cellchat.D9MOG.net[-which(cellchat.D9MOG.net$interaction_ST %in% cellchat.D9OVA.net$interaction_ST ),]


#all interactions in Day 7 OVA Scaf that are absent in Day 7 OVA Scaf
common_ccl_OVAD9 = cellchat.D9OVA.net.ccl.cxcl[which(cellchat.D9OVA.net.ccl.cxcl$interaction_ST %in% common_interactionsLate),]


#all CCL interactions that are the same Day 6 SC and Day 7 MOG Scaf
common_ccl_ScD10 = df10.net.ccl.cxcl[which(df10.net.ccl.cxcl$interaction_ST %in% common_interactionsLate),]


# CCL intereaactions in MOG and NOT OVA Day 9 scaf, and in Day 10 SC
common_ccl_ScD10_mog_notOVA = intersect(Mog_notinOVAD9$interaction_ST, common_ccl_ScD10$interaction_ST)
common_ccl_ScD10_mog_notOVA_net = cellchat.D9MOG.net.ccl.cxcl[which(cellchat.D9MOG.net.ccl.cxcl$interaction_ST %in% common_ccl_ScD10_mog_notOVA ),]


net.upMvO9$interaction_ST = paste0(net.upMvO9$interaction_name, "_",
                                   net.upMvO9$source, "_",
                                   net.upMvO9$target)

net.downMvO9$interaction_ST = paste0(net.downMvO9$interaction_name, "_",
                                     net.downMvO9$source, "_",
                                     net.downMvO9$target)

#CCL interactions in Day 10 SC and upregualted in diseased MOG day 9 Scaf (compared to Day 9 OVA)
common_ccl_ScD10_upmogvOVA = df10.net.ccl.cxcl[which( df10.net.ccl.cxcl$interaction_ST %in% net.upMvO9$interaction_ST),]

intersect(common_ccl_ScD10_upmogvOVA$interaction_ST, common_ccl_ScD10_mog_notOVA_net$interaction_ST)

#Figure shows CCL interactions in Day 10 SC and upregualted in diseased MOG day 9 Scaf (compared to Day 9 OVA)
png(file = "CommonCCL_Genechorddiagram_Scaf_late_upMogvOVA_inSC_colorcorrect.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "CommonCCL_Genechorddiagram_Scaf_late_upMogvOVA_inSC_colorcorrect.pdf",width = 8, height = 7)

netVisual_chord_gene(cellchat10,# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = "CCL",legend.pos.x = 8,
                     slot.name = 'net',net = common_ccl_ScD10_upmogvOVA,
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways in\n both the spinal cord and scaffold" ),
                     lab.cex = 1)
dev.off()

# CCL Pathways in MOG day 9 not in Day 10 SC
Dif_ccl_MOGD9 = cellchat.D9MOG.net.ccl.cxcl[-which(cellchat.D9MOG.net.ccl.cxcl$interaction_ST %in% common_ccl_ScD10_upmogvOVA$interaction_ST),]

# CCL Pathways in MOG day 9 not in Day 10 SC and absent in Day 9 OVA Scaf
Dif_ccl_MOGD9_notOVA = Dif_ccl_MOGD9[-which(Dif_ccl_MOGD9$interaction_ST %in% cellchat.D9OVA.net.ccl.cxcl$interaction_ST),]

# CCL Pathways in MOG day 9 not in Day 10 SC and not upregulated in Day 9 OVA Scaf
#Dif_ccl_MOGD9_notupOVA = Dif_ccl_MOGD9[-which(Dif_ccl_MOGD9$interaction_ST %in% net.downMvO9$interaction_ST),]
#Dif_ccl_MOGD9_notupOVA = Dif_ccl_MOGD9

# CCL Pathways in OVA and not in MOG scaf day 9 or  SC day 10
Dif_ccl_OVAD9 = cellchat.D9OVA.net.ccl.cxcl[-which(cellchat.D9OVA.net.ccl.cxcl$interaction_ST %in% common_interactionsLate),]

#CCL pathways in SC day 10 not in MOG scaf day 9
Dif_ccl_ScD10 = df10.net.ccl.cxcl[-which(df10.net.ccl.cxcl$interaction_ST %in% cellchat.D9MOG.net.ccl.cxcl$interaction_ST),]
#CCL pathways in SC day 10 not in MOG scaf day 9 or OVA scaf day 9
Dif_ccl_ScD101 = Dif_ccl_ScD10[-which(Dif_ccl_ScD10$interaction_ST %in% cellchat.D9OVA.net.ccl.cxcl$interaction_ST),]

#CCL pathways in SC day 10 not up in MOG scaf day 9 or down in OVA scaf day 9
Dif_ccl_ScD102 = df10.net.ccl.cxcl[-which(df10.net.ccl.cxcl$interaction_ST %in% net.upMvO9$interaction_ST),]
#Dif_ccl_ScD102 = Dif_ccl_ScD62[-which(Dif_ccl_ScD102$interaction_ST %in% net.downMvO9$interaction_ST),]

length(Dif_ccl_ScD62$interaction_ST)
length(Dif_ccl_ScD61$interaction_ST)

#figure shows CCL pathways in SC day 10 not in MOG scaf day 9 or OVA scaf day 9
pathways.show = "CCL"
png(file = "DifCCL_Genechorddiagram_SC_late_notinOVAoMOG_colorcorrect.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "DifCCL_Genechorddiagram_SC_late_notinOVAoMOG_colorcorrect.pdf",width = 8, height = 7)

netVisual_chord_gene(cellchat10,# sources.use = c(6) , targets.use = c(Tcells), 
                     #color.use = c("#8F539B", "#D03732", "#D43F89", "#4B7DB4", "#2B3270", "#E089B3", "blue","#9C5B32", "#E4993F", "darkgreen"),
                     signaling = "CCL",legend.pos.x = 8,
                     slot.name = 'net', net = Dif_ccl_ScD101, 
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways only in the spinal cord" ),
                     lab.cex = 1)
dev.off()


# CCL Pathways in MOG day 9 not in Day 10 SC and not upregulated in Day 9 OVA Scaf
png(file = "DifCCL_Genechorddiagram_Scaf_late_notinOVA_colorcorrect.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "DifCCL_Genechorddiagram_Scaf_late_notinOVA_colorcorrect.pdf",width = 8, height = 7)


netVisual_chord_gene(cellchat.D9MOG,# sources.use = c(6) , targets.use = c(Tcells), 
                     color.use = c("#377eb8", "#4daf4a", "#e4211c", "#bc9dcc", "grey", "#f781bf", "blue","#9C5B32", "#E4993F", "#54b0e4",
                                   "#b2df8a","#b2df8a","#a65629","#e3be00","darkgreen"),
                     signaling = c("CCL"),legend.pos.x = 8,
                     slot.name = 'net',net = Dif_ccl_MOGD9_notOVA, 
                     small.gap = 3, big.gap = 1, scale = T,
                     #sources.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     #targets.use = c("Neut","Mon","CD4 T","Mac","CD8 T","B Cells","DCs","NKs"),
                     title.name = paste("CCL","pathways only in the Scaffold" ),
                     lab.cex = 1)
dev.off()



#comparison of Scaf over time and SC over time ----

# Day 6 to Day 10 pos datset = day 10
net.up$interaction_ST = paste0(net.up$interaction_name, "_",
                               net.up$source, "_",
                               net.up$target)

net.down$interaction_ST = paste0(net.down$interaction_name, "_",
                               net.down$source, "_",
                               net.down$target)

net.upM7v9$interaction_ST = paste0(net.upM7v9$interaction_name, "_",
                                   net.upM7v9$source, "_",
                                   net.upM7v9$target)


net.downM7v9$interaction_ST = paste0(net.downM7v9$interaction_name, "_",
                                     net.downM7v9$source, "_",
                                     net.downM7v9$target)


commonup_ovrTime = intersect(net.up$interaction_ST, net.upM7v9$interaction_ST)

commondown_ovrTime = intersect(net.down$interaction_ST, net.downM7v9$interaction_ST)

