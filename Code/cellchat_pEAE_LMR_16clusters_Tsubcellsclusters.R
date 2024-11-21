

library(Seurat,lib.loc = .libPaths()[2])
library(SeuratObject,lib.loc = .libPaths()[2])
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

setwd("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/")

robj_path = "/Users/lailarad/Documents/ProgEAE_scRNAseq/LMR_16renamedclusters_wTcellsubcluster.rds"
save_path = "/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/"

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


data = RenameIdents(data, 
                    'Monocytes' = 'Mon',
                    "Macrophage" = "Mac",
                    "Complement Macrophage" = "Comp Mac",
                    "Neutrophil" = "Neut",
                    "NK Cell" = "NK",
                    "Stromal Cell" = "Stromal",
                    "CD11b+ DC" = "CD11b+ DC",
                    "Helper NK" = "Helper NK",
                    "B Cell" = "B Cell",
                    "CD103+ DC" = "CD103+ DC",
                    "Chemoattractant Monocyte" = "Inf Mon",
                    "Il4i1+ DC" = "Il4i1+ DC",
                    "Plasmablast" = "Plasmablast",
                    "Naive CD4" = "Naive CD4",
                    "CD4+ Effector Memory" = "CD4 Tem",
                    "Naive CD8" = "Naive CD8",
                    "Cytotoxic CD8 Tem" = "CD8 Tem",
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

plotPath = "/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/"

# D7 MOG Cellchat--------------------------------------------------
# Create CellChat Object
# input is a Seurat object
## Option 1: use the default cell identities of Seurat object

# sample: D7 MOG, D7 OVA, D9 MOG, D9 OVA
# Day: D7, D9
# Antigen: MOG, OVA 

Comparison = 'D7_MOG_'
Version = 'v1_'
Dataset = 'EAE_shortname_'

Idents(mog) = mog@meta.data$short
mog = SetIdent(mog, value = mog@meta.data$short)
cellchat <- createCellChat(object = mog, group.by = 'short', assay = "RNA")
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



# Preprocessing Expression Data for Cell-Cell Communication Analysis --------
# subset the expression data of signaling genes for saving computation cost

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
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# 
df.net.path <- subsetCommunication(cellchat, slot.name = "netP") # look at level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#cell.labels[c(2,3,4,17,18,19,20,21)] #T cells
#cell.labels[c(12,16)]

# df.indiv.T.source <- subsetCommunication(cellchat, sources.use = c(2,3,4,17,18,19,20,21) )#, targets.use = c(8,15) ) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.indiv.T.source, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Source_Individual', '.csv'), row.names = T)
# 
# df.indiv.T.rec <- subsetCommunication(cellchat, targets.use = c(3,9,10)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.indiv.T.rec, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Receiver_Individual', '.csv'), row.names = T)
# 
# df.net.T.source <- subsetCommunication(cellchat, sources.use =  c(3,9,10), slot.name = "netP") #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.net.T.source, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Source_Network', '.csv'), row.names = T)
# 
# df.net.T.rec <- subsetCommunication(cellchat, targets.use =  c(3,9,10), slot.name = "netP") #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.net.T.rec, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Receiver_Network', '.csv'), row.names = T)

# set to be the levels for T cell important communications 

# pathways = unique(df.net.T.rec$pathway_name, df.net.T.source$pathway_name)
# df.net <- subsetCommunication(cellchat, signaling = pathways) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.


# Infer cell-cell communication at signaling pathway level ----------------
cellchat <- computeCommunProbPathway(cellchat)


# Calculate aggregated communication network ------------------------------

cellchat <- aggregateNet(cellchat)

# visualize aggreggate network
groupSize <- as.numeric(table(cellchat@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

png(file = "D7MOGCommunication.png",units = "in",width = 6, height =6 , res = 400)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

# groupSize <- as.numeric(table(cellchat@idents))
# groupSize
# pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'interaction_weights', '.pdf'), onefile = TRUE)
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# dev.off()


mat <- cellchat@net$weight
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(5,4), mar = c(1, 1, 1, 1))
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

Tcells = c(14:19)

Bcells = c(9,13)


vertex.receiver = c(Tcells, Bcells)
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }

png(file = paste(plotPath, "AllpathwayscontributionLR_D7MOG.png",sep=""),
    units = "in",width = 30, height = 30, res = 400)
netAnalysis_contribution(cellchat, signaling = pathways.show.all,font.size = 5,
                         return.data = T)
dev.off()

netVisual_bubble(cellchat, sources.use = vertex.receiver, remove.isolate = FALSE)
netVisual_bubble(cellchat, targets.use = vertex.receiver, remove.isolate = FALSE)


# Part IV: Systems Analysis of cell-cell communcation network  ---------------------


# Compute and visualize the network centrality scores  --------------------
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[4], width = 8, height = 2.5, font.size = 10)
dev.off()

# Visualize the dominant senders (sources) and receivers (targets) --------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MHC-I", "MHC-II"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Identify signals contributing most to outgoing or incoming signa --------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 14)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-I", "MHC-II"))
ht

# (SLOW) Identify global communication patterns to explore how multiple c --------

library(NMF)
#> Loading required package: pkgmaker
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster
#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
#>   To enable shared memory capabilities, try: install.extras('
#> NMF
#> ')
#> 
#> Attaching package: 'NMF'
#> The following objects are masked from 'package:igraph':
#> 
#>     algorithm, compare
library(ggalluvial)

# v1_EAE_shortname_D7_MOG_outgoing_cluster_selection
k1 = selectK(cellchat, pattern = "outgoing")
k1
# both cophenetic and silhouette values begin to drop suddenly when number of outgoing patterns is 4 


# Comparison = 'Collective'

nPatterns = 4
dev.off()
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication', '.pdf'), onefile = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, height = 14)
dev.off()
# river plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication_riverplot', '.pdf'), onefile = TRUE)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication_dotplot', '.pdf'), onefile = TRUE)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

# Identify and visualize incoming communication pattern of target  --------
# v1_EAE_shortname_D7_MOG_incoming_cluster_selection
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_selection', '.pdf'), onefile = TRUE)
k2 = selectK(cellchat, pattern = "incoming")
k2
dev.off()

# could select 5 or 6 patterns, it's unclear 

nPatterns = 4
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_heatmap', '.pdf'), onefile = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, height = 15, width = 8)
dev.off()

# river plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_riverplot', '.pdf'), onefile = TRUE)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_dotplot', '.pdf'), onefile = TRUE)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

# Manifold and classification learning analysis of signaling netwo --------

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
#> 
# Enable parallelization
library(future)
plan()
plan("multisession", workers = 4)
trace(netClustering, edit=TRUE) #change multiprocess to multisession
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


# Identify signaling groups based on structure similarity -----------------
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


#chord----
table(cellchat@meta$ident)
df.net.T <- subsetCommunication(cellchat, sources.use = c(Tcells[1:4],6), targets.use = Tcells[1:4]) 
df.netP.T <- subsetCommunication(cellchat, sources.use = c(Tcells[1:4],6), targets.use = Tcells[1:4], slot.name = "netP") 
ggg2 = netVisual_chord_gene(cellchat, sources.use = Tcells,targets.use = Tcells, 
                            big.gap = 1 ,slot.name = 'net', net = df.net, lab.cex = 0.5, 
                            small.gap = 1, title.name = paste0("signaling in MOG Day 7"))

ggg2 = netVisual_chord_gene(cellchat, sources.use =  c(Tcells[1:4],6),targets.use =  c(Tcells[1:4],6), 
                            big.gap = 1 ,slot.name = 'netP', net = df.netP.T, lab.cex = 0.5, 
                            small.gap = 1, title.name = paste0("signaling in MOG Day 7"))
ggg2

netVisual_chord_gene(cellchat, sources.use =  c(Tcells[7], Bcells,6),targets.use =  c(Tcells[7], Bcells), 
                            big.gap = 4 ,slot.name = 'netP', net = df.net.path, lab.cex = 1, 
                            small.gap = 2, title.name = paste0("signaling in MOG Day 7"))

levels(cellchat@idents)
png(file = "D7MOGPathwaysTcellsAPCchorddiagram.png",units = "in",width = 6, height =6 , res = 400)
netVisual_chord_gene(cellchat, sources.use =  c(1,2,3,11),targets.use =  c(Tcells[1:4]), 
                     big.gap = 4 ,slot.name = 'netP', net = df.net.path, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("signaling in MOG Day 7"))

netVisual_chord_gene(cellchat, sources.use =  c(7,10,12, Tcells[1:4]),targets.use =  c(7,10,12, Tcells[1:4]), 
                     big.gap = 4 ,slot.name = 'netP', net = df.net.path, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("signaling in MOG Day 7"))

png(file = "D7MOGMHCchorddiagram.png",units = "in",width = 10, height = 6, res = 400)
netVisual_chord_gene(cellchat, sources.use =  c(1,2,3,11, Tcells),targets.use =  c(1,2,3,11,Tcells),
                     signaling = c("MHC-II", "MHC-I"),
                     big.gap = 2 ,slot.name = 'netP', net = df.net.path, lab.cex = 0.5, 
                     small.gap = 1, title.name = paste0("signaling in MOG Day 7"))

netVisual_chord_gene(cellchat,sources.use =  c(1,2,3,11, Tcells),targets.use =  c(1,2,3,11,Tcells),
                     signaling = c("MHC-II", "MHC-I"),
                     big.gap = 2 ,slot.name = 'netP', net = df.net.path, lab.cex = 0.5, 
                     small.gap = 1, title.name = paste0("signaling in MOG Day 7"))

netVisual_chord_gene(cellchat, sources.use =  c(1,2,3,11, Tcells[1:4]),targets.use =  c(1,2,3,11,Tcells[1:4]), 
                     signaling = c("MHC-I", "MHC-II"),
                     big.gap = 4 ,slot.name = 'net', net = df.net, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("MHC signaling in MOG Day 7"))

png(file = "D7MOGCCLCXCLchorddiagram.png",units = "in",width = 10, height = 6, res = 400)
netVisual_chord_gene(cellchat, sources.use =  c(1,2,3,11, Tcells),targets.use =  c(1,2,3,11,Tcells), 
                     signaling = c("CCL", "CXCL"),
                     big.gap = 4 ,slot.name = 'net', net = df.net, lab.cex = 1, 
                     small.gap = 2, title.name = paste0("CCL and CXCL signaling in MOG Day 7"))
dev.off()

netVisual_chord_gene(cellchat, #sources.use =  c(1,2,3,11),targets.use =  c(Tcells), 
                     signaling = c("TGFb"),
                     big.gap = 4 ,slot.name = 'net', net = df.net, lab.cex = 1, 
                     small.gap = 2, title.name = paste0("signaling in MOG Day 7"))
dev.off()
netVisual_chord_gene(cellchat,# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = c("MHC-II"),
                     slot.name = 'net',net = df.net, small.gap = 1, big.gap = 2,
                     title.name = paste(pathways.show, names(object.list)[i]),
                     lab.cex = 1)

netVisual_aggregate(cellchat, signaling = c("MHC-II"), layout = "chord", 
                    signaling.name = "MHC-II Pathway at Day 7 in EAE",
                    vertex.label.cex = 0.3,point.size = 0.5)

png(file = "D7MOGNK_TcellPathschorddiagram.png",units = "in",width = 10, height = 6, res = 400)
netVisual_chord_gene(cellchat, sources.use = c(5,8) , targets.use = c(15,17), 
                     #signaling = c("ITGAL-ITGB2"),
                     slot.name = 'netP',net = df.net.path, small.gap = 1, big.gap = 2,
                     title.name = "NK-Tem Signaling",
                     lab.cex = 0.8)

dev.off()

# Part V: Save the CellChat object ----------------------------------------
getwd()
unique(cellchat@meta$sample)



# D7 OVA Cellchat--------------------------------------------------
# Create CellChat Object
# input is a Seurat object
## Option 1: use the default cell identities of Seurat object

# sample: D7 MOG, D7 OVA, D9 MOG, D9 OVA
# Day: D7, D9
# Antigen: MOG, OVA 

Comparison = 'D7_OVA_'
Version = 'v1_'
Dataset = 'EAE_shortname_'

Idents(ova) = ova@meta.data$short
ova = SetIdent(ova, value = ova@meta.data$short)
cellchat <- createCellChat(object = ova, group.by = 'short', assay = "RNA")
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



# Preprocessing Expression Data for Cell-Cell Communication Analysis --------
# subset the expression data of signaling genes for saving computation cost

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
#cellchat <- projectData(cellchat, PPI.mouse)

# Compute Communicaiton probability and infer network ---------------------
cellchat <- computeCommunProb(cellchat)
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 5)

# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# 
df.net.path <- subsetCommunication(cellchat, slot.name = "netP") # look at level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.

#cell.labels[c(2,3,4,17,18,19,20,21)] #T cells
#cell.labels[c(12,16)]

# df.indiv.T.source <- subsetCommunication(cellchat, sources.use = c(3,9,10) )#, targets.use = c(8,15) ) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.indiv.T.source, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Source_Individual', '.csv'), row.names = T)
# 
# df.indiv.T.rec <- subsetCommunication(cellchat, targets.use = c(3,9,10)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.indiv.T.rec, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Receiver_Individual', '.csv'), row.names = T)
# 
# df.net.T.source <- subsetCommunication(cellchat, sources.use =  c(3,9,10), slot.name = "netP") #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.net.T.source, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Source_Network', '.csv'), row.names = T)
# 
# df.net.T.rec <- subsetCommunication(cellchat, targets.use =  c(3,9,10), slot.name = "netP") #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.net.T.rec, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Receiver_Network', '.csv'), row.names = T)
# 
# # set to be the levels for T cell important communications 
# 
# pathways = unique(df.net.T.rec$pathway_name, df.net.T.source$pathway_name)
# df.net <- subsetCommunication(cellchat, signaling = pathways) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.


# Infer cell-cell communication at signaling pathway level ----------------
cellchat <- computeCommunProbPathway(cellchat)


# Calculate aggregated communication network ------------------------------

cellchat <- aggregateNet(cellchat)

# visualize aggreggate network
groupSize <- as.numeric(table(cellchat@idents))
groupSize
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


groupSize <- as.numeric(table(cellchat@idents))
groupSize
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'interaction_weights', '.pdf'), onefile = TRUE)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


mat <- cellchat@net$weight
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(5,5), mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

cellchat@netP$pathways

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
cellchat@netP$MOG$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)

Tcells = c(14:19)

Bcells = c(9,13)

vertex.receiver = c(Tcells)
for (i in 1:length(pathways.show.all)) {
  i = 5
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  #ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

png(file = paste(plotPath, "AllpathwayscontributionLR.png",sep=""),
    units = "in",width = 30, height = 30, res = 400)
netAnalysis_contribution(cellchat, signaling = pathways.show.all, 
                         font.size = 5,
                         return.data = T)
dev.off()

netVisual_bubble(cellchat, sources.use = c(Tcells, Bcells), 
                 remove.isolate = FALSE, font.size = 5)
netVisual_bubble(cellchat, targets.use = c(Tcells, Bcells), remove.isolate = FALSE)


# Part IV: Systems Analysis of cell-cell communcation network  ---------------------


# Compute and visualize the network centrality scores  --------------------
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[3], width = 8, height = 2.5, font.size = 10)
dev.off()

# Visualize the dominant senders (sources) and receivers (targets) --------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("MHC-I", "MHC-II"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Identify signals contributing most to outgoing or incoming signa --------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", 
                                         width = 10, height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",
                                         width = 10, height = 20)
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-I", "MHC-II"))
ht

# (SLOW) Identify global communication patterns to explore how multiple c --------

library(NMF)
#> Loading required package: pkgmaker
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster
#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
#>   To enable shared memory capabilities, try: install.extras('
#> NMF
#> ')
#> 
#> Attaching package: 'NMF'
#> The following objects are masked from 'package:igraph':
#> 
#>     algorithm, compare
library(ggalluvial)

# v1_EAE_shortname_D7_MOG_outgoing_cluster_selection
k1 = selectK(cellchat, pattern = "outgoing")
k1
# both cophenetic and silhouette values begin to drop suddenly when number of outgoing patterns is 4 


# Comparison = 'Collective'

nPatterns = 4
dev.off()
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication', '.pdf'), onefile = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing",
                                          k = nPatterns,
                                          height = 15)
dev.off()
# river plot
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication_riverplot', '.pdf'), onefile = TRUE)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication_dotplot', '.pdf'), onefile = TRUE)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

# Identify and visualize incoming communication pattern of target  --------
# v1_EAE_shortname_D7_MOG_incoming_cluster_selection
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_selection', '.pdf'), onefile = TRUE)
k2 = selectK(cellchat, pattern = "incoming")
k2
dev.off()



nPatterns = 4

#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_heatmap', '.pdf'), onefile = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, height = 15)
dev.off()

# river plot
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_riverplot', '.pdf'), onefile = TRUE)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
#pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_dotplot', '.pdf'), onefile = TRUE)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

# Manifold and classification learning analysis of signaling netwo --------

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
#> 
# Enable parallelization
library(future)
plan()
plan("multisession", workers = 4)
trace(netClustering, edit=TRUE) #change multiprocess to multisession
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


# Identify signaling groups based on structure similarity -----------------
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

#chord----
table(cellchat@meta$ident)
df.net.T <- subsetCommunication(cellchat, sources.use = c(Tcells[1:4],6), targets.use = Tcells[1:4]) 
df.netP.T <- subsetCommunication(cellchat, sources.use = c(Tcells[1:4],6), targets.use = Tcells[1:4], slot.name = "netP") 
ggg2 = netVisual_chord_gene(cellchat, sources.use = Tcells,targets.use = Tcells, 
                            big.gap = 1 ,slot.name = 'net', net = df.net, lab.cex = 0.5, 
                            small.gap = 1, title.name = paste0("signaling in MOG Day 7"))

ggg2 = netVisual_chord_gene(cellchat, sources.use =  c(Tcells[1:4],6),targets.use =  c(Tcells[1:4],6), 
                            big.gap = 1 ,slot.name = 'netP', net = df.netP.T, lab.cex = 0.5, 
                            small.gap = 1, title.name = paste0("signaling in MOG Day 7"))
ggg2

netVisual_chord_gene(cellchat, sources.use =  c(Tcells[4], Bcells,6),targets.use =  c(Tcells[4], Bcells), 
                     big.gap = 4 ,slot.name = 'netP', net = df.net.path, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("signaling in OVA Day 7"))


# Part V: Save the CellChat object ----------------------------------------
getwd()
unique(cellchat@meta$sample)
saveRDS(cellchat, file = "cellchat_D7OVA.rds")

# Comparison of D7 MOG vs D7 OVA -----
#comparison analysis ----
setwd("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/D7 MvO/")
cellchat.D7MOG <- readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/cellchat_D7MOG.rds")
cellchat.D7OVA <- readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/cellchat_D7OVA.rds")
cellchat.D7MOG = updateCellChat(cellchat.D7MOG)
cellchat.D7OVA = updateCellChat(cellchat.D7OVA)

object.list <- list(OVA = cellchat.D7OVA, MOG = cellchat.D7MOG)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
cellchat@meta
unique(cellchat@meta$sample)
object.list$OVA@meta$ident
object.list$MOG@meta$ident
unique(cellchat@meta$ident)
Tcells

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2) )
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",  )
gg1 + gg2
dev.off()

#Differential number of interactions or interaction strength among different cell populations
#where red (or blue) colored edges represent increased (or decreased) signaling 
#in the second dataset MOG compared to the first one OVA.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)

#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


pdf(file = "NumInteractionsD7.pdf", 
    width = 16, height = 7, pointsize = 12)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

dev.off()

length(cell.labels)
length()
unique(object.list$OVA@meta$short)
unique(object.list$MOG@meta$short)
unique(object.list$OVA@meta$ident)
unique(object.list$MOG@meta$ident)

group.cellType <- c(rep("B Cell", 3), rep("CD8 T", 4),  rep("Il4i1+ DC", 4), rep("Plasmablast",4))
group.cellType <- factor(group.cellType, levels = c("B Cell","CD8 T","Il4i1+ DC","Plasmablast"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD8 Tem", signaling.exclude = "MIF", top.label = 0.25)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4 Tem", signaling.exclude = c("MIF"), top.label = 0.25)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Plasmablast", signaling.exclude = c("MIF"), top.label = 0.25)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "B Cell", signaling.exclude = c("MIF"), top.label = 0.25)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2, gg3, gg4))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = c("functional"), label.size = 3.5, pathway.remove.show = F)
#> 2D visualization of signaling networks from datasets 1 2

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

rankSimilarity(cellchat, type = "functional")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.10.0
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> The new InteractiveComplexHeatmap package can directly export static 
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

Tcells = c(14:19)

Bcells = c(9,13)


levels(cellchat@meta$short)

netVisual_bubble(cellchat, signaling = c("CD52", "CD23", "TGFb"), comparison = c(1, 2))
pairLR <- extractEnrichedLR(cellchat, signaling = "CD23", geneLR.return = FALSE)
pairLR
netVisual(cellchat, pairLR.use= pairLR[1] , vertex.receiver= c(Tcells) )
netVisual(cellchat, signaling = pathway.union[5], layout = "hierarchy")


netVisual_bubble(cellchat,targets.use = c(Tcells, Bcells), sources.use = c(Tcells, Bcells), comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(3,8,9,10,15),targets.use = c(3,8,9,10,15), comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, targets.use = c(3,9,10), comparison = c(1, 2), angle.x = 45)
allcomp = netVisual_bubble(cellchat,  comparison = c(1, 2), angle.x = 45)
table = allcomp$data
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat, targets.use = c(Tcells), sources.use = c(1,2), comparison = c(1, 2), font.size.title = 24, max.dataset = 2, title.name = "Increased signaling in MOG", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, targets.use = c(Tcells), sources.use = c(1,2),  comparison = c(1, 2), font.size.title = 24, max.dataset = 1, title.name = "Decreased signaling in MOG", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, targets.use = c(3,9,10), comparison = c(1, 2), font.size.title = 24, max.dataset = 2, title.name = "Increased signaling in MOG", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, targets.use =  c(3,9,10),  comparison = c(1, 2), font.size.title = 24, max.dataset = 1, title.name = "Decreased signaling in MOG", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2


# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MOG"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in MOG
net.up <- subsetCommunication(cellchat, net = net, datasets = "MOG",ligand.logFC = 0.2)
#net.down <- subsetCommunication(cellchat, net = net, datasets = "OVA",ligand.logFC = 0.2)


#write.csv(net.up,paste(getwd(),"/net.up_Day7_MogvsOVA.csv", sep = ""), row.names = T)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in OVA, i.e.,downregulated in MOG
net.down <- subsetCommunication(cellchat, net = net, datasets = "OVA",ligand.logFC = -0.2, receptor.logFC = -0.2)

#write.csv(net.down,paste(getwd(),"/net.down_Day7_MogvsOVA.csv", sep = ""), row.names = T)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
c(3,8,9,10,15)
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, targets.use = c(Tcells, Bcells), 
                        sources.use = c(Tcells, Bcells),  comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg1
netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(3,8,9,10,15), targets.use = c(3,8,9,10,15),
                 comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
ggg1 = netVisual_chord_gene(object.list[[2]], sources.use = c(3,8,9,10,15),targets.use = c(3,8,9,10,15), 
                     big.gap = 1 ,slot.name = 'net', net = net.up, lab.cex = 0.5, 
                     small.gap = 1, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[2]], sources.use = c(2,4),targets.use = c(3,9,10), 
                     big.gap = 1 ,slot.name = 'net', net = net.up, lab.cex = 0.5, 
                     small.gap = 1, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

ggg2 = netVisual_chord_gene(object.list[[2]], sources.use = c(3,8,9,10,15),targets.use = c(3,8,9,10,15), 
                     big.gap = 1 ,slot.name = 'net', net = net.down, lab.cex = 0.5, 
                     small.gap = 1, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "CD86", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, features = "Il2rb", split.by = "datasets", colors.ggplot = T) #CD122-CD8+ T cells pathogenic
plotGeneExpression(cellchat, features = "Il17f", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, features = "Rorc", split.by = "datasets", colors.ggplot = T) 
plotGeneExpression(cellchat, features = "Ifng", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, features = "Fcer1a", split.by = "datasets", colors.ggplot = T)


#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(3,9,10),  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(3,8,9,10,15), targets.use = c(3,8,9,10,15),
                 comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

#> Comparing communications on a merged object
gg1 + gg2
ggg1 + ggg2


gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, targets.use = c(3,9,10), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, targets.use = c(3,9,10),  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

computeEnrichmentScore(net.down, species = 'mouse')
computeEnrichmentScore(net.up, species = 'mouse')

# Chord diagram error
par(mfrow = c(1,2), xpd=TRUE, mar = c(7,7,7,7))
object.list[[2]]@meta$short

png(file = "chorddiagram.png",units = "in",width = 8, height = 6, res = 400)
netVisual_chord_gene(object.list[[2]] , targets.use = c(7,10,13), sources.use = Tcells,
                     big.gap = 3 ,slot.name = 'net', net = net.up, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

png(file = "chorddiagram.png",units = "in",width = 8, height = 6, res = 400)
netVisual_chord_gene(object.list[[2]], targets.use = c(3,9,10), 
                     big.gap = 3 ,slot.name = 'net', net = net.up, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

netVisual_chord_gene(object.list[[1]], sources.use = Tcells, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

png(file = "TregMHCIIchorddiagram.png",units = "in",width = 10, height = 6, res = 400)

pathways.show <- c("MHC-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

MCHIID7MVO= netVisual_bubble(cellchat, signaling = pathways.show, 
                 comparison = c(1,2), angle.x = 90,
                 line.on = T, remove.isolate = F,
                 title.name = "MHC-II D7 MOG vs. OVA",
                 return.data = T
)

view(MCHIID7MVO$communication)


par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("MHC-I") 
pathways.show <- c("CCL")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pathways.show <- c("CD52") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CD23") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
c(IL1, ITGAL-ITGB2, CD86,CD80, CXCL, CCL, PD-L1,CSF)
png(file = "IL1_Il18chorddiagram.png",units = "in",width = 10, height = 6, res = 400)
pathways.show <- c("IL1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


c(IL1, ITGAL-ITGB2, CD86,CD80, CXCL, CCL, PD-L1,CSF)
png(file = "IL1_Il18chorddiagram.png",units = "in",width = 10, height = 6, res = 400)
pathways.show <- c("CD86") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

table(cellchat@meta$ident)

pathways.show <- c("MHC-II") 
pathways.show <- c("CD23") 
pathways.show <- c("IL1") 

#pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
#i =2
#png(file = "IL1_Il18Genechorddiagram.png",units = "in",width = 10, height = 6, res = 400)
netVisual_chord_gene(object.list[[i]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net, small.gap = 1, big.gap = 2,
                     title.name = paste(pathways.show, names(object.list)[i]),
                     lab.cex = 1)
}
dev.off()

pathways.show <- c("CD22")
#pathways.show <- unique(net.up.Treg$pathway_name)
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), 
                           attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, unique(object.list$MOG@meta$Day) ,names(object.list)[i]))
}

cellchat@meta

pathways.show <- c("THBS", "FN1 ")
i =2

pathways.show.all <- c(cellchat@netP$MOG$pathways, cellchat@netP$OVA$pathways)

pdf(file ="cellchat.pdf", width = 25, height =16)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  #i =2
  # png(file = "IL1_Il18Genechorddiagram.png",units = "in",width = 10, height = 6, res = 400)
  netVisual_chord_gene(object.list[[i]], sources.use = c("Treg") , targets.use = c("CD4 Tem",
                                                                                   "CD8 Tem",
                                                                                   "Naive CD4",
                                                                                   "Naive CD8",
                                                                                   "gd", "Treg"), 
                       signaling = pathways.show.all,legend.pos.x = 8,
                       slot.name = 'net',net = net, small.gap = 2, big.gap = 2,
                       title.name = paste(names(object.list)[i]),
                       lab.cex = 1)
}
dev.off()
getwd()

weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show.all) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show.all, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

cellchat@meta$short
levels(cellchat@idents$joint)
net.up.Treg = net.up[which(net.up$source == "Treg" | net.up$target == "Treg"),]
net.down.Treg = net.down[which(net.down$source == "Treg" | net.down$target == "Treg"),]

netVisual_chord_gene(object.list[[2]] ,
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'net',net = net.down.Treg, small.gap = 3, big.gap = 1,
                     title.name = paste("Treg pathways downregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

gene.up.Treg <- extractGeneSubsetFromPair(net.up.Treg, cellchat)
gene.down.Treg <- extractGeneSubsetFromPair(net.down.Treg, cellchat)

pairLR.use.up.Treg = net.up.Treg[, "interaction_name", drop = F]

netVisual_bubble(cellchat, pairLR.use = pairLR.use.up.Treg, 
                        comparison = c(1, 2),  angle.x = 90, 
                        remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

net.up.Treg.ECM_R = net.up.Treg[which(net.up.Treg$annotation == "ECM-Receptor"),]
net.down.Treg.ECM_R = net.down.Treg[which(net.down.Treg$annotation == "ECM-Receptor"),]



netVisual_chord_gene(object.list[[2]] ,
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'net',net = net.up.Treg.ECM_R, small.gap = 3, big.gap = 1,
                     title.name = paste("Treg ECM-Receptor pathways upregulated in" ,names(object.list)[2], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

netVisual_chord_gene(object.list[[2]] , 
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'net',net = net.down.Treg.ECM_R, small.gap = 3, big.gap = 1,
                     title.name = paste("Treg ECM-Receptor pathways downregulated in" ,names(object.list)[2], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

netVisual_chord_gene(object.list[[2]] ,
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'netP',
                     net = net.up[which(net.up$annotation == "ECM-Receptor"),], 
                     small.gap = 3, big.gap = 1,
                     title.name = paste("ECM-Receptor pathways upregulated in" ,names(object.list)[2], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

netVisual_chord_gene(object.list[[2]] ,
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'netP',
                     net = net.down[which(net.down$annotation == "ECM-Receptor"),], 
                     small.gap = 3, big.gap = 1,
                     title.name = paste("ECM-Receptor pathways downregulated in" ,names(object.list)[2], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

dev.off()

netVisual_chord_gene(object.list[[2]] , 
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'net',
                     net = net.down.Treg.ECM_R, 
                     small.gap = 3, big.gap = 1,
                     title.name = paste("pathways downregulated in" ,names(object.list)[2], unique(cellchat@meta$Day) ),
                     lab.cex = 1)



netVisual_chord_gene(object.list[[2]] ,
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'netP',net = net.up.Treg, small.gap = 3, big.gap = 1,
                     title.name = paste("Treg pathways upregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

netVisual_chord_gene(object.list[[2]], targets.use = c(19), 
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'netP',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste("pathways upregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

netVisual_chord_gene(object.list[[2]] , sources.use = c(19), 
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'net',net = net.down, small.gap = 3, big.gap = 1,
                     title.name = paste("pathways downregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

netVisual_chord_gene(object.list[[2]], sources.use = c(19), 
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste("pathways upregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

pathways.show <- c("CCL") 

pairLR.use.up.CCL = net.up[which(net.up$pathway_name == "CCL"), "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up.CCL, #sources.use = 4, targets.use = c(5:11), 
                        comparison = c(1,2),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Day 7: Up-regulated signaling in ", names(object.list)[2]))

pdf(file = "CCL_Dotplot_Upreg_MogvsOVA_Day7.pdf", 
    width = 25, height = 5, pointsize = 12)
gg1

dev.off()
getwd()
png(file = "CCL_Genechorddiagram_Upreg_MogvsOVA_Day7.png",units = "in",width = 8, height = 7, res = 400)

library(svglite)
svglite::svglite(file = "CCL_Genechorddiagram_Upreg_MogvsOVA_Day7_svg.svg", 
                 width = 8, height = 7, pointsize = 12)


pdf(file = "CCL_Genechorddiagram_Upreg_MogvsOVA_Day7.pdf", 
    width = 8, height = 7, pointsize = 12)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show, legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated in" ,names(object.list)[2], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

dev.off()

unique(cellchat@meta$ident)

Tcells = c(6, 8, 9, 10, 15, 18)
netVisual_chord_gene(object.list[[1]], sources.use = c(18) , targets.use = c(Tcells), 
                     signaling = pathways.show.all,legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

setwd("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/D7 MvO/")

png(file = "CCL_Genechorddiagram_Downreg_MogvsOVA_Day7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.down, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

png(file = "CCL_Genechorddiagram_Downreg_MogvsOVA_Day7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)

netVisual_chord_gene(object.list[[i]], sources.use = c(5,8,15,17) , targets.use = c(5,8,15,17), 
                    legend.pos.x = 8,
                     slot.name = 'net',net = net, small.gap = 1, big.gap = 2,
                     title.name = paste(pathways.show, names(object.list)[i]),
                     lab.cex = 0.5)


png(file = "NKTcellGenechorddiagram.png",units = "in",width = 10, height = 6, res = 400)
netVisual_chord_gene(object.list[[i]], sources.use = c(5,8) , targets.use = c(15,17), 
                     legend.pos.x = 8, signaling = pathway.union,
                     slot.name = 'netP',net = net.up, small.gap = 1, big.gap = 2,
                     title.name = paste(pathways.show, names(object.list)[i]),
                     lab.cex = 0.5)

netVisual_chord_gene(object.list[[2]] , targets.use = c(Tcells), sources.use = c(6),
                     big.gap = 3 ,slot.name = 'net', net = net, lab.cex = 0.5, 
                     
                     small.gap = 2, title.name = paste0("Up-regulated signaling in ", 
                                                        names(object.list)[2]))


pathways.show <- c("PARs") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

 
# Chord diagram
pathways.show <- c("CD52") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "chord", 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      vertex.label.cex = 0.3,point.size = 0.5)
}
dev.off()

netVisual_chord_gene(cellchat, sources.use =  c(Bcells),#targets.use =  c(Tcells), 
                     big.gap = 1 ,slot.name = 'net', net = net.up, lab.cex = 0.1, 
                     small.gap = 1, title.name = paste0("signaling in MOG Day 7"))


getwd()
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_D7MOGvOVA.rds")


# D9 MOG Cellchat--------------------------------------------------
# Create CellChat Object
# input is a Seurat object
## Option 1: use the default cell identities of Seurat object

# sample: D7 MOG, D7 OVA, D9 MOG, D9 OVA
# Day: D7, D9
# Antigen: MOG, OVA 
setwd("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/")
D9 = subset(data, idents = 'D9')

Idents(D9) = D9@meta.data$Antigen
mog = subset(D9, idents = "MOG")
ova = subset(D9, idents = "OVA")


Comparison = 'D9_MOG_'
Version = 'v1_'
Dataset = 'EAE_shortname_'

Idents(mog) = mog@meta.data$short
mog = SetIdent(mog, value = mog@meta.data$short)
cellchat <- createCellChat(object = mog, group.by = 'short', assay = "RNA")
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



# Preprocessing Expression Data for Cell-Cell Communication Analysis --------
# subset the expression data of signaling genes for saving computation cost

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
#cellchat <- projectData(cellchat, PPI.mouse)

# Compute Communicaiton probability and infer network ---------------------
cellchat <- computeCommunProb(cellchat)
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# 
df.net.path <- subsetCommunication(cellchat, slot.name = "netP") # look at level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
levels(cellchat@idents)
Tcells = c(14:19)
Bcells = c(9,13)

# df.indiv.T.source <- subsetCommunication(cellchat, sources.use = c(3,9,10) )#, targets.use = c(8,15) ) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.indiv.T.source, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Source_Individual', '.csv'), row.names = T)
# 
# df.indiv.T.rec <- subsetCommunication(cellchat, targets.use = c(3,9,10)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.indiv.T.rec, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Receiver_Individual', '.csv'), row.names = T)
# 
# df.net.T.source <- subsetCommunication(cellchat, sources.use =  c(3,9,10), slot.name = "netP") #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.net.T.source, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Source_Network', '.csv'), row.names = T)
# 
# df.net.T.rec <- subsetCommunication(cellchat, targets.use =  c(3,9,10), slot.name = "netP") #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.net.T.rec, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Receiver_Network', '.csv'), row.names = T)
# 
# # set to be the levels for T cell important communications 
# 
# pathways = unique(df.net.T.rec$pathway_name, df.net.T.source$pathway_name)
# df.net <- subsetCommunication(cellchat, signaling = pathways) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.


# Infer cell-cell communication at signaling pathway level ----------------
cellchat <- computeCommunProbPathway(cellchat)


# Calculate aggregated communication network ------------------------------

cellchat <- aggregateNet(cellchat)

# visualize aggreggate network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


groupSize <- as.numeric(table(cellchat@idents))
groupSize
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'interaction_weights', '.pdf'), onefile = TRUE)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


mat <- cellchat@net$weight
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(5,4), mar = c(1, 1, 1, 1))
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

vertex.receiver = Tcells
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }

png(file = paste(plotPath, "AllpathwayscontributionLR.png",sep=""),
    units = "in",width = 30, height = 30, res = 400)
netAnalysis_contribution(cellchat, signaling = pathways.show.all)
dev.off()

netVisual_bubble(cellchat, sources.use = Tcells, remove.isolate = FALSE)
netVisual_bubble(cellchat, targets.use = Tcells, remove.isolate = FALSE)


# Part IV: Systems Analysis of cell-cell communcation network  ---------------------


# Compute and visualize the network centrality scores  --------------------
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[3], width = 8, height = 2.5, font.size = 10)
dev.off()

# Visualize the dominant senders (sources) and receivers (targets) --------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Identify signals contributing most to outgoing or incoming signa --------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("MHC-I", "MHC-II"))
ht

# (SLOW) Identify global communication patterns to explore how multiple c --------

library(NMF)
#> Loading required package: pkgmaker
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster
#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
#>   To enable shared memory capabilities, try: install.extras('
#> NMF
#> ')
#> 
#> Attaching package: 'NMF'
#> The following objects are masked from 'package:igraph':
#> 
#>     algorithm, compare
library(ggalluvial)

# v1_EAE_shortname_D7_MOG_outgoing_cluster_selection
k = selectK(cellchat, pattern = "outgoing")
k
# both cophenetic and silhouette values begin to drop suddenly when number of outgoing patterns is 4 


# Comparison = 'Collective'

nPatterns = 3
dev.off()
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication', '.pdf'), onefile = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, height = 12)
dev.off()
# river plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication_riverplot', '.pdf'), onefile = TRUE)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication_dotplot', '.pdf'), onefile = TRUE)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

# Identify and visualize incoming communication pattern of target  --------
# v1_EAE_shortname_D7_MOG_incoming_cluster_selection
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_selection', '.pdf'), onefile = TRUE)
k1 = selectK(cellchat, pattern = "incoming")
k1
dev.off()

# could select 5 or 6 patterns, it's unclear 

nPatterns = 4
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_heatmap', '.pdf'), onefile = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, height = 15)
dev.off()

# river plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_riverplot', '.pdf'), onefile = TRUE)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_dotplot', '.pdf'), onefile = TRUE)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

# Manifold and classification learning analysis of signaling netwo --------

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
#> 
# Enable parallelization
library(future)
plan()
plan("multisession", workers = 4)
trace(netClustering, edit=TRUE) #change multiprocess to multisession
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


# Identify signaling groups based on structure similarity -----------------
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

#chord----
levels(cellchat@meta$ident)
df.net.T <- subsetCommunication(cellchat, sources.use = c(Tcells[1:4],6), targets.use = Tcells[1:4]) 
df.netP.T <- subsetCommunication(cellchat, sources.use = c(Tcells[1:4],6), targets.use = Tcells[1:4], slot.name = "netP") 
ggg2 = netVisual_chord_gene(cellchat, sources.use = Tcells,targets.use = Tcells, 
                            big.gap = 1 ,slot.name = 'net', net = df.net, lab.cex = 0.5, 
                            small.gap = 1, title.name = paste0("signaling in MOG Day 9"))

ggg2 = netVisual_chord_gene(cellchat, sources.use =  c(Tcells[1:4],6),targets.use =  c(Tcells[1:4],6), 
                            big.gap = 1 ,slot.name = 'netP', net = df.netP.T, lab.cex = 0.5, 
                            small.gap = 1, title.name = paste0("signaling in MOG Day 9"))
ggg2

netVisual_chord_gene(cellchat, sources.use =  c(Tcells[7], Bcells,6),targets.use =  c(Tcells[7], Bcells), 
                     big.gap = 4 ,slot.name = 'netP', net = df.net.path, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("signaling in MOG Day 9"))



# Part V: Save the CellChat object ----------------------------------------
getwd()
cellchat@meta$sample
saveRDS(cellchat, file = "cellchat_D9MOG_LS.rds")

# D9 OVA Cellchat--------------------------------------------------
# Create CellChat Object
# input is a Seurat object
## Option 1: use the default cell identities of Seurat object

# sample: D7 MOG, D7 OVA, D9 MOG, D9 OVA
# Day: D7, D9
# Antigen: MOG, OVA 


setwd("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/13clusters/D9 MvO")
Comparison = 'D9_OVA_'
Version = 'v1_'
Dataset = 'EAE_shortname_'

Idents(ova) = ova@meta.data$short
ova = SetIdent(ova, value = ova@meta.data$short)
cellchat <- createCellChat(object = ova, group.by = 'short', assay = "RNA")
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



# Preprocessing Expression Data for Cell-Cell Communication Analysis --------
# subset the expression data of signaling genes for saving computation cost

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
#cellchat <- projectData(cellchat, PPI.mouse)

# Compute Communicaiton probability and infer network ---------------------
cellchat <- computeCommunProb(cellchat)
#> triMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2022-11-26 08:34:38]"
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-26 08:35:36]"
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 5)

# Extract inferred network as data frame ----------------------------------
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# 
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# 
df.net.path <- subsetCommunication(cellchat, slot.name = "netP") # look at level of signaling pathways

#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
levels(cellchat@meta$short)
levels(cellchat@idents)
Tcells = c(14:19)
Bcells = c(9,13)

# df.indiv.T.source <- subsetCommunication(cellchat, sources.use = c(3,9,10) )#, targets.use = c(8,15) ) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.indiv.T.source, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Source_Individual', '.csv'), row.names = T)
# 
# df.indiv.T.rec <- subsetCommunication(cellchat, targets.use = c(3,9,10)) #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.indiv.T.rec, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Receiver_Individual', '.csv'), row.names = T)
# 
# df.net.T.source <- subsetCommunication(cellchat, sources.use =  c(3,9,10), slot.name = "netP") #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.net.T.source, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Source_Network', '.csv'), row.names = T)
# 
# df.net.T.rec <- subsetCommunication(cellchat, targets.use =  c(3,9,10), slot.name = "netP") #gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# write.csv(df.net.T.rec, file = paste0(plotPath, Version, Dataset, Comparison, 'TCell_Receiver_Network', '.csv'), row.names = T)
# 
# # set to be the levels for T cell important communications 
# 
# pathways = unique(df.net.T.rec$pathway_name, df.net.T.source$pathway_name)
# df.net <- subsetCommunication(cellchat, signaling = pathways) #gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
# 

# Infer cell-cell communication at signaling pathway level ----------------
cellchat <- computeCommunProbPathway(cellchat)


# Calculate aggregated communication network ------------------------------

cellchat <- aggregateNet(cellchat)

# visualize aggreggate network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


groupSize <- as.numeric(table(cellchat@idents))
groupSize
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'interaction_weights', '.pdf'), onefile = TRUE)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


mat <- cellchat@net$weight
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'signaling_clusters', '.pdf'), onefile = TRUE)
par(mfrow = c(5,4), mar = c(1, 1, 1, 1))
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

# vertex.receiver = c(3,9,10)
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }

png(file = paste(plotPath, "AllpathwayscontributionLR.png",sep=""),
    units = "in",width = 30, height = 30, res = 400)
netAnalysis_contribution(cellchat, signaling = pathways.show.all)
dev.off()

netVisual_bubble(cellchat, sources.use = cell.labels[c(3,9,10)], remove.isolate = FALSE)
netVisual_bubble(cellchat, targets.use = cell.labels[c(3,9,10)], remove.isolate = FALSE)


# Part IV: Systems Analysis of cell-cell communcation network  ---------------------


# Compute and visualize the network centrality scores  --------------------
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[1], width = 8, height = 2.5, font.size = 10)
dev.off()

# Visualize the dominant senders (sources) and receivers (targets) --------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL","CD52"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2


# Identify signals contributing most to outgoing or incoming signa --------
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))


# (SLOW) Identify global communication patterns to explore how multiple c --------

library(NMF)
#> Loading required package: pkgmaker
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster
#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
#>   To enable shared memory capabilities, try: install.extras('
#> NMF
#> ')
#> 
#> Attaching package: 'NMF'
#> The following objects are masked from 'package:igraph':
#> 
#>     algorithm, compare
library(ggalluvial)

# v1_EAE_shortname_D7_MOG_outgoing_cluster_selection
k = selectK(cellchat, pattern = "outgoing")
k
# both cophenetic and silhouette values begin to drop suddenly when number of outgoing patterns is 4 


# Comparison = 'Collective'

nPatterns = 2
dev.off()
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication', '.pdf'), onefile = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()
# river plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication_riverplot', '.pdf'), onefile = TRUE)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication_dotplot', '.pdf'), onefile = TRUE)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

# Identify and visualize incoming communication pattern of target  --------
# v1_EAE_shortname_D7_MOG_incoming_cluster_selection
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_selection', '.pdf'), onefile = TRUE)
k1 = selectK(cellchat, pattern = "incoming")
k1
dev.off()

# could select 5 or 6 patterns, it's unclear 

nPatterns = 2
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_heatmap', '.pdf'), onefile = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()

# river plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_riverplot', '.pdf'), onefile = TRUE)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

#> Please make sure you have load `library(ggalluvial)` when running this function

# dot plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'incoming_communication_dotplot', '.pdf'), onefile = TRUE)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

# Manifold and classification learning analysis of signaling netwo --------

cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
#> 
# Enable parallelization
library(future)
plan()
plan("multisession", workers = 4)
trace(netClustering, edit=TRUE) #change multiprocess to multisession
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)

# netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)


# Identify signaling groups based on structure similarity -----------------
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural", umap.method = "uwot")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


# Part V: Save the CellChat object ----------------------------------------
getwd()
cellchat@meta$sample
saveRDS(cellchat, file = "cellchat_D9OVA_LS.rds")

# Comparison of D9 MOG vs D9 OVA -----
#comparison analysis ----
setwd("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/D9 MvO/")
cellchat.D9MOG <- readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/cellchat_D9MOG_LS.rds")
cellchat.D9OVA <- readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/cellchat_D9OVA_LS.rds")
cellchat.D9MOG = updateCellChat(cellchat.D9MOG)
cellchat.D9OVA = updateCellChat(cellchat.D9OVA)
object.list <- list(OVA = cellchat.D9OVA, MOG = cellchat.D9MOG)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat
cellchat@meta
unique(cellchat@meta$sample)
unique(cellchat@meta$ident)
unique(object.list$OVA@meta$ident)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Differential number of interactions or interaction strength among different cell populations
#where red (or blue) colored edges represent increased (or decreased) signaling 
#in the second dataset compared to the first one.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

pdf(file = "NumInteractionsD9MvO.pdf", 
    width = 16, height = 7, pointsize = 12)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

dev.off()

length(cell.labels)

levels(unique(object.list$OVA@meta$short))

group.cellType <- c(rep("Imm DC", 3), rep("DC", 3), rep("CD4 EM", 3), rep("Cytotoxic CD8", 3), rep("Naive CD4", 3), rep("CD8 Tcm",3),rep("Treg",3))
group.cellType <- factor(group.cellType, levels = c("Imm DC", "DC", "CD4 EM", "Cytotoxic CD8", "Naive CD4","CD8 Tcm", "Treg"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) +
    ylim(0,15) + xlim(0,20)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Plasma", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "DC", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4 EM", top.label = 1)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cytotoxic CD8", top.label = 1)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg5 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Th17 T", top.label = 1)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2, gg3, gg4))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

rankSimilarity(cellchat, type = "functional")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

pathways.up.mog.list = c("CD137","GALECTIN", "FGF", "CD45", "CADM",
                         "SELL", "CSF", "CD22", "PERIOSTIN", "ANGPTL",
                         "NECTIN", "HSPG", "APRIL", "PECAM1", "FASLG",
                         "CHEMERIN", "GAS", "TGFb")

library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.10.0
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> The new InteractiveComplexHeatmap package can directly export static 
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

stim_mol = c( "MHC-I", "MHC-II","CD80", "CD86")
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = stim_mol, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling =  stim_mol, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

netVisual_bubble(cellchat, sources.use = c(3,9,10), comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, targets.use = c(3,9,10), comparison = c(1, 2), angle.x = 45)
allcomp = netVisual_bubble(cellchat,  comparison = c(1, 2), angle.x = 45)
table = allcomp$data

netVisual_bubble(cellchat, signaling = stim_mol, comparison = c(1, 2), angle.x = 45)

#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat, sources.use = c(3,9,10), comparison = c(1, 2), font.size.title = 24, max.dataset = 2, title.name = "Increased signaling in MOG", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use =  c(3,9,10),  comparison = c(1, 2), font.size.title = 24, max.dataset = 1, title.name = "Decreased signaling in MOG", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, targets.use = c(3,9,10), comparison = c(1, 2), font.size.title = 24, max.dataset = 2, title.name = "Increased signaling in MOG", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, targets.use =  c(3,9,10),  comparison = c(1, 2), font.size.title = 24, max.dataset = 1, title.name = "Decreased signaling in MOG", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "MOG"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in MOG

net.up <- subsetCommunication(cellchat, net = net, datasets = "MOG",ligand.logFC = 0.2)
#write.csv(net.up,paste(getwd(),"/net.up_Day9_MogvsOVA.csv", sep = ""), row.names = T)


# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in OVA, i.e.,downregulated in MOG
net.down <- subsetCommunication(cellchat, net = net, datasets = "OVA",ligand.logFC = -0.1, receptor.logFC = -0.1)
#write.csv(net.down,paste(getwd(),"/net.down_Day9_MogvsOVA.csv", sep = ""), row.names = T)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)



pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(3,9,10), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(3,9,10),  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, targets.use = c(3,9,10), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, targets.use = c(3,9,10),  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

computeEnrichmentScore(net.down, species = 'mouse')
computeEnrichmentScore(net.up, species = 'mouse')

# Chord diagram error
par(mfrow = c(1,2), xpd=TRUE, mar = c(7,7,7,7))
object.list[[2]]@meta$short

png(file = "chorddiagram.png",units = "in",width = 8, height = 6, res = 400)
netVisual_chord_gene(object.list[[2]], targets.use = Tcells, 
                     big.gap = 3 ,slot.name = 'net', net = net.up, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

netVisual_chord_gene(object.list[[1]], targets.use = Tcells, 
                     big.gap = 3 ,slot.name = 'net', net = net.down, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))


netVisual_chord_gene(object.list[[2]], sources.use = DCs, 
                     big.gap = 3 ,slot.name = 'net', net = net.up, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

netVisual_chord_gene(object.list[[2]], sources.use = c(Tcells,Bcell), targets.use =c(Tcells,Bcell), 
                     big.gap = 3 ,slot.name = 'net', net = net.up, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(Tcells,Bcell), targets.use =c(Tcells,Bcell), 
                     big.gap = 3 ,slot.name = 'net', net = net.down, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))



netVisual_chord_gene(object.list[[1]],  targets.use =Tcells, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

pathways.show <- c("CD22")
pathways.show <- c("FASLG")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("MHC-I") 
pathways.show <- c("SELL")
pathways.show <- c("CD23")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

# Chord diagram from D7----
pathways.show <- c("MHC-I") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pathways.show <- c("CD52") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("CD23") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
c(IL1, ITGAL-ITGB2, CD86,CD80, CXCL, CCL, PD-L1,CSF)
#png(file = "IL1_Il18chorddiagram.png",units = "in",width = 10, height = 6, res = 400)
pathways.up.mog.list
pathways.show <- c( "LAMININ") 
pathways.show <- c( "NECTIN") 
pathways.show <- c( "ICAM") 
pathways.show <- c( "ITGAL-ITGB2") 
pathways.show <- c( "MHC-II") 

pdf(file = "ICAMD9MvO.pdf", 
    width = 16, height = 7, pointsize = 12)
pathways.show <- c( "ICAM") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i], unique(cellchat@meta$Day)))
}

dev.off()

pdf(file = "ICAMD9MvO_dotplot.pdf", 
    width = 16, height = 7, pointsize = 12)

netVisual_bubble(cellchat, signaling = pathways.show, 
                 comparison = c(2,1), angle.x = 90,
                 line.on = T, remove.isolate = T,
                 #max.dataset = 1, 
                 title.name = "ICAMD9MvO",
                 return.data = F
                )

sub_net.up = net.up[which(net.up$pathway_name == pathways.show),]
pairLR.use.up = sub_net.up[, "interaction_name", drop = F]

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        comparison = c(1,2),  angle.x = 90, 
                        remove.isolate = F,
                        title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
sub_net.down = net.down[which(net.down$pathway_name == pathways.show),]
pairLR.use.down = sub_net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,
                        comparison = c(1,2),  angle.x = 90, 
                        remove.isolate = F,
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2


ggsave(file = "ICAMD9MvO_dotplot.svg", width = 16, height = 7, pointsize = 12)
dev.off()

pdf(file = "ITGAL_ITGB2D9MvO.pdf", 
    width = 16, height = 7, pointsize = 12)
pathways.show <- c( "ITGAL-ITGB2") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i], unique(cellchat@meta$Day)))
}

dev.off()

netVisual_bubble(cellchat, signaling = pathways.show, 
                 comparison = c(2,1), angle.x = 90,
                 line.on = T, remove.isolate = T,
                 title.name = "ITGAL_ITGB2D9MvO",
                 return.data = F
)

ggsave(file = "ITGAL_ITGB2D9MvO_dotplot.svg", width = 16, height = 7, pointsize = 12)

pdf(file = "LamininD9MvO.pdf", 
    width = 16, height = 7, pointsize = 12)
pathways.show <- c( "LAMININ")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i], unique(cellchat@meta$Day)))
}

dev.off()

pairLR.use.up = net.up[, "interaction_name", drop = F]
netVisual_bubble(cellchat, signaling = pathways.show, 
                 comparison = c(1,2), angle.x = 90,
                 line.on = T, remove.isolate = F,
                 title.name = "LamininD9MvO",
                 return.data = T
)

ggsave(file = "LamininD9MvO_dotplot.svg", width = 16, height = 7, pointsize = 12)


#genes of interest in inflam th1 cells Cd226, Cd96, Cd44
#cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "CD226", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, features = "Cd44", split.by = "datasets", colors.ggplot = T) #CD122-CD8+ T cells pathogenic
plotGeneExpression(cellchat, features = "Cd226", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, features = "Cd96", split.by = "datasets", colors.ggplot = T) 
plotGeneExpression(cellchat, features = "Nectin1", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, features = "Bdkrb2", split.by = "sample", colors.ggplot = T)
plotGeneExpression(cellchat, features = "Il10", split.by = "sample", colors.ggplot = T)

stim_mol = c("MHC-I", "MHC-II","CD80","CD86")
plotGeneExpression(cellchat, signaling = stim_mol, split.by = "sample", colors.ggplot = T)


table(cellchat@meta$ident)

pathways.show <- c("MHC-II") 
pathways.show <- c("CD23") 
pathways.show <- c("IL1") 
#pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  #i =2
 # png(file = "IL1_Il18Genechorddiagram.png",units = "in",width = 10, height = 6, res = 400)
  netVisual_chord_gene(object.list[[i]],# sources.use = c(6) , targets.use = c(Tcells), 
                       signaling = pathways.show,legend.pos.x = 8,
                       slot.name = 'net',net = net, small.gap = 1, big.gap = 2,
                       title.name = paste(pathways.show, names(object.list)[i]),
                       lab.cex = 1)
}
dev.off()



netVisual_chord_gene(object.list[[i]], sources.use = c(5,8,15,17) , targets.use = c(5,8,15,17), 
                     legend.pos.x = 8,
                     slot.name = 'net',net = net, small.gap = 1, big.gap = 2,
                     title.name = paste(pathways.show, names(object.list)[i]),
                     lab.cex = 0.5)
#png(file = "NKTcellGenechorddiagram.png",units = "in",width = 10, height = 6, res = 400)
levels(cellchat@meta$ident)
netVisual_chord_gene(object.list[[i]], sources.use = c(12) , targets.use = c(15,17), 
                     signaling = pathway.union[20],
                     slot.name = 'net',net = net, small.gap = 1, big.gap = 2,
                     title.name = paste(pathways.show, names(object.list)[i]),
                     lab.cex = 0.5)

netVisual_chord_gene(object.list[[i]], sources.use = c(Tcells, Bcells), targets.use = c(Tcells, Bcells), 
                     signaling = pathway.union,
                     slot.name = 'netP',net = net, small.gap = 1, big.gap = 2,
                     #title.name = paste(pathways.show, names(object.list)[i]),
                     lab.cex = 0.5)

netVisual_chord_gene(object.list[[2]] ,targets.use = c(Tcells, Bcells),
                     big.gap = 3 ,slot.name = 'net', net = net.up, lab.cex = 0.5, 
                     signaling = c("CCL"),
                     small.gap = 2, title.name = paste0("Up-regulated signaling in ", 
                                                        names(object.list)[2]))
dev.off()
getwd()
pathways.show <- c("CCL") 
png(file = "CCL_Genechorddiagram_Upreg_MogvsOVA_Day9.png",units = "in",width = 8, height = 7, res = 400)

pdf(file = "CCL_Genechorddiagram_Upreg_MogvsOVA_Day9.pdf", 
    width = 8, height = 7, pointsize = 12)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)
dev.off()

pairLR.use.up.CCL = net.up[which(net.up$pathway_name == "CCL"), "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up.CCL, #sources.use = 4, targets.use = c(5:11), 
                        comparison = c(1,2),  angle.x = 90, remove.isolate = T,
                        title.name = paste0("Day 9: Up-regulated signaling in ", names(object.list)[2]))

pdf(file = "CCL_Dotplot_Upreg_MogvsOVA_Day9.pdf", 
    width = 25, height = 5, pointsize = 12)
gg1

dev.off()

pathways.show <- c("CD86") 
png(file = "CCL_Genechorddiagram_Downreg_MogvsOVA_Day9.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.down, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list)[i], unique(cellchat@meta$Day) ),
                     lab.cex = 1)
dev.off()

pathways.show <- c("PARs") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


# Chord diagram
pathways.show <- c("CD52") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "chord", 
                      signaling.name = paste(pathways.show, names(object.list)[i]),
                      vertex.label.cex = 0.3,point.size = 0.5)
}
dev.off()

pathways.show <- stim_mol
j = 2
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = stim_mol[j], layout = "chord", 
                      signaling.name = paste(stim_mol[j], names(object.list)[i]),
                      vertex.label.cex = 0.3,point.size = 0.5)
}

netVisual_chord_gene(cellchat, sources.use =  c(Bcells),#targets.use =  c(Tcells), 
                     big.gap = 1 ,slot.name = 'net', net = net.up, lab.cex = 0.1, 
                     small.gap = 1, title.name = paste0("signaling in MOG Day 7"))

dev.off()
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_D9MOGvOVA.rds")



# Comparison of D9 MOG vs D7 MOG -----
#comparison analysis ----

setwd("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/D7vD9 Mog/")
cellchat.D9MOG <- readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/cellchat_D9MOG_LS.rds")
cellchat.D7MOG <- readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/16clusters_plusTcells_BMES/cellchat_D7MOG.rds")
cellchat.D9MOG = updateCellChat(cellchat.D9MOG)
cellchat.D7MOG = updateCellChat(cellchat.D7MOG)
object.list <- list(D7 = cellchat.D7MOG, D9 = cellchat.D9MOG)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#cellchat <- readRDS("/Users/lailarad/Documents/ProgEAE_scRNAseq/CellChat/cellchat_comparisonAnalysis_MOGD9vD7.rds")
cellchat
cellchat@meta
unique(cellchat@meta$sample)
unique(cellchat@meta$short)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

#Differential number of interactions or interaction strength among different cell populations
#where red (or blue) colored edges represent increased (or decreased) signaling 
#in the second dataset compared to the first one.

#mistake in source code fixed: https://github.com/sqjin/CellChat/issues/604
#paste in line 109: igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
trace(netVisual_diffInteraction, edit=TRUE) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count")
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

length(cell.labels)

unique(object.list$D7@meta$short)

group.cellType <- c(rep("Imm DC", 3), rep("DC", 3), rep("CD4 EM", 3), rep("Cytotoxic CD8", 3), rep("Naive CD4", 3), rep("CD8 Tcm",3),rep("Treg",3))
group.cellType <- factor(group.cellType, levels = c("Imm DC", "DC", "CD4 EM", "Cytotoxic CD8", "Naive CD4","CD8 Tcm", "Treg"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

group.cellType <- c(rep("Imm DC", 2), rep("DC", 2), rep("CD4 T", 2), rep("CD8 T", 2), rep("Th17 T", 2), rep("NK",2), rep("Fib", 2), rep("Neut",1), rep("Mac", 1))
group.cellType <- factor(group.cellType, levels = c("Imm DC", "DC", "CD4 T", "CD8 T", "Th17 T", "NK", "Fib", "Neut","Mac"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "count.merged", label.edge = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight.merged", label.edge = T)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) +
    ylim(0,20) + xlim(0,20)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Plasmablast", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "DC", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg3 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Naive CD4", top.label = 1)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
gg4 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cytotoxic CD8", top.label = 1)
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2,gg3, gg4))
dev.off()

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

rankSimilarity(cellchat, type = "functional")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

library(ComplexHeatmap)
#> Loading required package: grid
#> ========================================
#> ComplexHeatmap version 2.10.0
#> Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
#> Github page: https://github.com/jokergoo/ComplexHeatmap
#> Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
#> 
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
#>   genomic data. Bioinformatics 2016.
#> 
#> The new InteractiveComplexHeatmap package can directly export static 
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#> 
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

cell.labels[c(3,9,10)] #T cells
cell.labels[c(8,15)]

netVisual_bubble(cellchat,comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat, sources.use = c(3,9,10),comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in D9", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use =  c(3,9,10),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in D9", angle.x = 45, 
                        remove.isolate = T, line.on = T)
#> Comparing communications on a merged object
gg1 + gg2

group.cellType <- c(rep("Imm DC", 2), rep("DC", 2), rep("CD4 T", 2), rep("CD8 T", 2), rep("Th17 T", 2), rep("NK",2), rep("Fib", 2), rep("Neut",1), rep("Mac", 1))
group.cellType <- factor(group.cellType, levels = c("Imm DC", "DC", "CD4 T", "CD8 T", "Th17 T", "NK", "Fib", "Neut","Mac"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

netVisual_bubble(cellchat, sources.use = c(3,9,10), comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, targets.use = c(3,9,10), comparison = c(1, 2), angle.x = 45)
allcomp = netVisual_bubble(cellchat,  comparison = c(1, 2), angle.x = 45)
table = allcomp$data
#> Comparing communications on a merged object

gg1 <- netVisual_bubble(cellchat, sources.use = c(3,9,10), comparison = c(1, 2), font.size.title = 24, max.dataset = 2, title.name = "Increased signaling in D9", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use =  c(3,9,10),  comparison = c(1, 2), font.size.title = 24, max.dataset = 1, title.name = "Decreased signaling in D9", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, targets.use = c(3,9,10), comparison = c(1, 2), font.size.title = 24, max.dataset = 2, title.name = "Increased signaling in D9", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, targets.use =  c(3,9,10),  comparison = c(1, 2), font.size.title = 24, max.dataset = 1, title.name = "Decreased signaling in D9", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "D9"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in D9
net.up <- subsetCommunication(cellchat, net = net, datasets = "D9",ligand.logFC = 0.2)
write.csv(net.up,paste(plotPath,"/net.up_Day9v7_Mog.csv", sep = ""), row.names = T)

# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in D7, i.e.,downregulated in D9
net.down <- subsetCommunication(cellchat, net = net, datasets = "D7",ligand.logFC = -0.1, receptor.logFC = -0.1)
write.csv(net.down,paste(plotPath,"/net.down_Day9v7_Mog.csv", sep = ""), row.names = T)


gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(3,9,10), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(3,9,10),  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, targets.use = c(3,9,10), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, targets.use = c(3,9,10),  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

computeEnrichmentScore(net.down, species = 'mouse')
computeEnrichmentScore(net.up, species = 'mouse')

# Chord diagram error
par(mfrow = c(1,2), xpd=TRUE, mar = c(7,7,7,7))
object.list[[2]]@meta$short

Tcells = c(14:19)
Bcells = c(9,13)

png(file = "chorddiagram.png",units = "in",width = 8, height = 6, res = 400)
netVisual_chord_gene(object.list[[2]], targets.use = Tcells , 
                     big.gap = 3 ,slot.name = 'net', net = net.up, lab.cex = 0.5, 
                     small.gap = 2, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
dev.off()

netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

pathways.show <- c("MHC-I") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
pathways.show <- c("MHC-I") 
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

dev.off()


---
computeEnrichmentScore(net.down, species = 'mouse')
computeEnrichmentScore(net.up, species = 'mouse')

# Chord diagram error
par(mfrow = c(1,2), xpd=TRUE, mar = c(1,1,1,1))
object.list[[2]]@meta$short

#object.list[[2]]= D9 MOG
netVisual_chord_gene(object.list[[2]], targets.use = c(10), big.gap = 2 ,slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], targets.use = c(3), slot.name = 'net', net = net.down, lab.cex = 0.8, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
dev.off()

pathways.show <- c("MHC-I") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# Chord diagram
pathways.show <- c("CCL") 
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pathways.show <- c("CD52") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

pathways.show <- c("MHC-I", "MHC-II", "CD80", "CD86") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))


# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# Chord diagram
#pathways.show <- c("CCL") 
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

nPatterns = 4
dev.off()
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication', '.pdf'), onefile = TRUE)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()
# river plot
pdf(file = paste0(plotPath, Version, Dataset, Comparison, 'outgoing_communication_riverplot', '.pdf'), onefile = TRUE)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()


pathway.union
pathways.show <- c("CCL") 
png(file = "CCL_Genechorddiagram_Upreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated in" ,names(object.list)[i] ),
                     lab.cex = 1)
dev.off()

png(file = "CCL_Genechorddiagram_Downreg_Mog_Day9v7.png",units = "in",width = 8, height = 7, res = 400)
netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.down, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways downregulated in" ,names(object.list)[i] ),
                     lab.cex = 1)
dev.off()

saveRDS(cellchat, file = "cellchat_comparisonAnalysis_MOGD9vD7.rds")

pathways.show <- c( "ICAM") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}


pathways.show <- c( "ITGAL-ITGB2") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}
getwd()
write.csv(net.up,"net.up.Day9v7.csv")
write.csv(net.down,"net.down.Day9v7.csv")


pathways.show <- c( "CCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                      edge.weight.max = weight.max[1], edge.width.max = 5, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     signaling = "CCL",legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated in" ,names(object.list)[2]),
                     lab.cex = 1)
dev.off()

pathways.show <- c( "CCL") 
getwd()
png(file = "CCL_Genechorddiagram_Up_Scaf_overtime.png",units = "in",width = 8, height = 7, res = 400)
pdf(file = "CCL_Genechorddiagram_Up_Scaf_overtime.pdf",width = 8, height = 7)

netVisual_chord_gene(object.list[[2]],# sources.use = c(6) , targets.use = c(Tcells), 
                     color.use = c("#D03732","#4B7DB4", "#2B3270","#8F539B", "#E4993F","blue","#9C5B32",
                                               "#E4993F","blue","#5B3B1B","#B31E1A",
                                   "blue","blue","blue",
                                   "#D43F89", "blue","blue",
                                               "#E089B3", "#D43F89",
                                   "blue"),
                     signaling = pathways.show,legend.pos.x = 8,
                     slot.name = 'net',net = net.up, small.gap = 3, big.gap = 1,
                     title.name = paste(pathways.show,"pathways upregulated \nin the scaffold over time" ),
                     lab.cex = 1)

#ggsave(file = "CCL_Genechorddiagram_Up_Scaf_overtime.svg", width = 8, height = 7, pointsize = 12)

dev.off()

