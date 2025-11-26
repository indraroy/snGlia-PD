#this is the template I used for the kamath dataset
#replace the objct name for each separte obj
library("ggplot2")
library("Seurat")
library("cowplot")
library("clustree")
library(patchwork)
library(dplyr)
library("Matrix")
library("MAST")
# read in obj
kamath_mg <- readRDS("kamath_mg.Rds")
colnames(kamath_mg@meta.data)

#make sure you have the right assay & set ident
DefaultAssay(kamath_mg) <- 'RNA'
Idents(kamath_mg) <- kamath_mg@meta.data$disease__ontology_label

#DGE analysis
kamath_mg_DGE <- FindMarkers(kamath_mg, test.use = "MAST", ident.1 = "Parkinson disease", ident.2 = "normal")
write.csv(kamath_mg_DGE, "kamath_mg PD markers.csv", row.names = TRUE)
#for pretty volcano plot I used a 0 logfc threshold
kamath_mg_DGE <- FindMarkers(kamath_mg, test.use = "MAST", ident.1 = "Parkinson's disease", ident.2 = "normal", logfc.threshold = 0)

##GEORGINA'S CODE FOR VOLCANO PLOTS

#gseGO analysis
library(BiocManager)
library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#bring in saved DGE list
kamath_mgDGElist <- read.csv("kamath_mg PD markers.csv")
#prep GO pipeline
kamath_mg_genelist <- kamath_mgDGElist$avg_log2FC
names(kamath_mg_genelist) <- rownames(kamath_mgDGElist)
kamath_mg_genelist<-na.omit(kamath_mg_genelist)
kamath_mg_genelist = sort(kamath_mg_genelist, decreasing = TRUE)
## gene enrichment time
kamath_mgGOenrichment <- gseGO(geneList=kamath_mg_genelist, 
             ont ="MF", 
             keyType = "SYMBOL",
             minGSSize = 2, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05,
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db)

# save as csv results
mg_GOres.df <- as.data.frame(kamath_mgGOenrichment@result)
write.csv(mg_GOres.df, "kamath_mg_gseGO_MF results.csv")

## GEORGINA'S CODE FOR GO PLOTS