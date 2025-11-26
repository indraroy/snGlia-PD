#################################################################################
#                                                                               #
# R code for plots on paper:                                                    #
# Roy I, Thomas RA, Jimenez Ambriz G, Farhan S, Piscopo VEC, Durcan TM, 2025.   #
# Single nuclear RNA sequencing shows altered microglial and astrocytic         #
# functions in post-mortem Parkinson’s disease tissue. (In Revision)            #
#                                                                               #
#################################################################################

#############################
## VOLCANO PLOTS
#############################

library(dplyr)
library(tidyverse)
library(DOSE)
library(org.Hs.eg.db)
library(ggrepel)

# 1. Genes to highlight
lab_genes_Kamath_mg=c("HSP90AB1", "HSP90AA1", "HSPA6", "BAG3", "UTY", "GPNMB", "DNAJA4", "LRRK2",
                      "IL1B", "TIMP1", "TMEM163", "HLA-DRA", "CX3CR1", 
                      "HIF3A")

# 2. Load DEG list
DEG<-read.csv("Kamath_mg_DEG_list.csv")

# 3. Change the lower ceiling p value (p_val=0) to p_val = 3.521222e-315
DEG$log10_pval_adj<--log10(DEG$p_val_adj)
DEG$log10_pval_adj[DEG$log10_pval_adj==Inf]<-500

# 4. Add a column (diffexpressed) of NAs
DEG$diffexpressed <- "not significant"

# 5. Set diffexpressed column to "SIG" if log2Foldchange > 0.5 and pvalue < 0.05, or
# if log2Foldchange < -0.5 and pvalue < 0.05
DEG$diffexpressed[(DEG$avg_log2FC > 0.5 & DEG$p_val_adj < 0.05) |
                          (DEG$avg_log2FC < -0.5 & DEG$p_val_adj < 0.05)] <- "SIG"

# 6. Set diffexpressed column to "SIG" if the gene is in the markers list
DEG$diffexpressed[DEG$X%in%lab_genes_Kamath_mg] <- "Markers"

# 7.Create a new column "delabel2" that will contain the name of genes differentially
# expressed (NA in case they are not).
DEG<-DEG %>%  mutate(delabel2 = case_when(X%in%lab_genes_Kamath_mg ~ X))

# 8. Plot
options(ggrepel.max.overlaps = Inf)
mycolors <- c("black","gray","red")
ggplot(data=DEG, aes(x=avg_log2FC, y=log10_pval_adj, col=diffexpressed, label=delabel2)) +
  geom_point(data = DEG[!is.na(DEG$delabel2),], color = "red") +
  geom_point(alpha = 0.8, size=2) + 
  theme_minimal() +
  geom_text_repel(point.size = 4, # data point size
                  size = 3, 
                  color="black",
                  point.padding = 0, # additional padding around each point
                  min.segment.length = 0, # draw all line segments
                  box.padding = 0.5) + # additional padding around each text label)
  geom_point(data = DEG[!is.na(DEG$delabel2),], color = "black") +
  scale_color_manual(values= mycolors) +
  ggtitle("Kamath microglia PD volcano")+
  geom_vline(xintercept=c(-0.5, 0.5), col="gray37", linetype="dotted") +
  geom_hline(yintercept=-log10(0.05), col="gray37", linetype=3)


#################################################
# GSEA gene set enrichment analysis
# DOTPLOTS AND CNETPLOTS
#################################################
library(tidyverse)
library(DOSE)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(fgsea)
library(data.table)
library(ggplot2)
library("cowplot")
library(dplyr)
library("Matrix")

# 1. Load DEG list 
DEG<-read.csv("Kamath_mg_DEG_list.csv")
head(DEG)

# 2. GO pipeline preparation

# 2.1 Extract folchange and genes
# Make a named vector of fold change values
DEG_genelist <- DEG$avg_log2FC
names(DEG_genelist) <- DEG$X
DEG_genelist<-na.omit(DEG_genelist)
DEG_genelist = sort(DEG_genelist, decreasing = TRUE)
head(DEG_genelist)
length(DEG_genelist)

# 2.2 Load gene annotations
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)


# 3. Run GO Gene Set Enrichment Analysis

gGOenrichment_MF <- gseGO(geneList=DEG_genelist,
                          ont ="MF",
                          keyType = "SYMBOL",
                          minGSSize = 2,
                          maxGSSize = 1000,
                          pvalueCutoff = 0.1,
                          verbose = TRUE,
                          OrgDb = org.Hs.eg.db,
                          eps=0)

# 4. Check up at number of terms
dim(gGOenrichment_MF)

# 5. Make Dotplot
dotplot(gGOenrichment_MF, showCategory=10, font.size=7, split=".sign",
        label_format=40) + facet_grid(.~.sign) + 
  ggtitle("kamath MG MF diff_functional analysis") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold"),
        legend.text = element_text(size=7), legend.title = element_text(size=8))


# 6. Make cnetplots

# 6.1 Pick first 5 up-regulated terms and first 5 down regulated terms 
GO_up <- head(gGOenrichment_MF$Description[gGOenrichment_MF$NES > 0], 5)
GO_down<-head(gGOenrichment_MF$Description[gGOenrichment_MF$NES < 0], 5)
GO_cnet_terms<-c(GO_up, GO_down)

# 6.2 Make Cnetplot
x2 <- pairwise_termsim(gGOenrichment_MF) 
cnetplot(x2,
         showCategory =GO_cnet_terms,
         node_label="all",
         layout =  'kk',
         color.params = list(foldChange = DEG_genelist),
         cex.params = list(category_node = 0.6, gene_node = 1, category_label = 0.5,
                           gene_label = 0.4),
         shadowtext ="category") + 
  theme(legend.text=element_text(size=8))+
  theme(legend.key.size = unit(0.5,"line"))


#######################
# HEATMAPS
#######################

library(tidyverse)
library('Seurat')
library(ComplexHeatmap)

##########################################################################
# PART 1: This part of the analysis was run in an interactive session on a 
# HPC Alliance Canada Cluster.

# 1. Read Seurat object
seu <- readRDS("~/runs/georgia/scRNA_downtream_analysis/indra_DGE/kamath_Dec-2024/kamath_mg.Rds")

# 2. Add a column to metadata to relabel disease ontology label
seu$disease_ontology_label<-seu$disease__ontology_label
seu$disease_ontology_label[seu$disease__ontology_label == 'normal'] <- 'control'
seu$disease_ontology_label[seu$disease__ontology_label == 'Parkinson disease'] <- "Parkinson's disease"

# 3. Make sure you have the right assay
DefaultAssay(seu) <- 'RNA'

# 4. Filter unwanted cells
Idents(seu) <- "Cell_Subtype"
cells_to_remove <- WhichCells(seu, idents = "Macro_CD200R1")
seu_filtered <- subset(seu, cells = cells_to_remove, invert = TRUE)

# 5. Set slot and idents
DefaultAssay(seu_filtered) <- 'RNA'
Idents(seu_filtered) <-seu_filtered@meta.data$disease_ontology_label

# 6. Scale data data
seu_filtered <- ScaleData(seu_filtered)

# 7. Check cell number
ncol(seu_filtered)
table(Idents(seu_filtered))

# 8. Remove excluded samples from analysis and downsample
seu_filtered_5kCells <- subset(seu_filtered, subset = disease_ontology_label !=
                                 "Lewy body dementia")
set.seed(111)
sampled.cells <- sample(x =colnames(seu_filtered_5kCells), size = 5000,
                        replace = F)
seu_filtered_5kCells <-  subset(seu_filtered_5kCells, cells = sampled.cells)

# 9. List all genes of interest
genes<-c("HSP90AB1", "BAG3", "DNAJA4", "HSPA1B", "NAMPT", "STIP1", "CHORDC1", "FKBP4",
         "GPNMB", "AIF1", "CTSD", "HIF3A", "C1QA", "CD83", "HLA-DRA", "CSF3R",
         "CSF2RA", "MEG3", "IL4R", "IL17RA", "LILRB4", "CD74","CX3CR1", "P2RY12",
         "TGFB1", "CD68", "TMEM163","LRRK2", "CLEC7A", "APOE",
         "LPL", "TREM2", "ITGAX", "CTSB")

# 10. Extract interesting genes from scaled gene expression matrix
seu_filtered_5kCells<-seu_filtered_5kCells$RNA@scale.data[mg_genes,]

# 11. Extract metadata from filtered object
annotations<-seu_filtered_5kCells[[]]

# 12. Write extracted matrices to csv file
write.table(seu_filtered_5kCells, file = "mat_kamath_mg_filtered_5kCells.csv",
            sep = ",", row.names = TRUE, col.names = TRUE)
write.table(annotation, file = "mat_kamath_mg_filtered_annotations.csv",
            sep = ",", row.names = TRUE, col.names = TRUE)



########################################################################
# PART 2: This part of the analysis was run in a personal computer.

# 1. Read filtered scaled expression matrix and annotations
gene_expression<-read.table(file = "mat_kamath_mg_filtered_5kCells.csv",
                                                   sep= ",", row.names = 1)
annotations<-read.table(file = "mat_kamath_mg_filtered_annotations.csv",
                                                  sep= ",", row.names = 1)

# 2. Check up that all cell names are in annotations
head(colnames(gene_expression), n=2)
head (rownames(annotations), n=2)
# Rename cell names in annotations to be the same than in expression matrix
rownames(annotations) <- gsub("-", ".", rownames(annotations))
all(colnames(gene_expression) %in% rownames(annotations))


# 3. Make vectors of gene sets function
Chaperone_protein_binding_functions<-c("HSP90AB1", "BAG3", "DNAJA4", "HSPA1B",
                                       "NAMPT", "STIP1", "CHORDC1", "FKBP4")

Immune_binding_functions<-c("GPNMB", "AIF1", "CTSD", "HIF3A", "C1QA", "CD83",
                            "HLA-DRA", "CSF3R","CSF2RA", "MEG3", "IL4R", "IL17RA",
                            "LILRB4", "CD74")

Microglial_state_markers<-c("CX3CR1", "P2RY12", "TGFB1", "CD68", "TMEM163")

Damage_disease_associated<-c("LRRK2", "CLEC7A", "APOE", "LPL", "TREM2", "ITGAX","CTSB")

genes_2<-c(Chaperone_protein_binding_functions,Immune_binding_functions,
           Microglial_state_markers,Damage_disease_associated)


# 4. Create data frame of genes function
Chaperone_pbf <- rep("Chaperone protein binding functions", 8)
Immune_bf <- rep("Immune binding functions", 14)
Microglial_sm <- rep("Microglial state markers", 5)
Damage_da <- rep("Damage disease associated", 7)
Function_mg <- c(Chaperone_pbf,Immune_bf,Microglial_sm,Damage_da)
names(Function_mg) <- genes_2
df_function <- as.data.frame(Function_mg)

# 5. Sort annotations and expression matrix by status and cell type
annotations<-annotations[order(annotations$disease_ontology_label,
                               annotations$Cell_Subtype), ]
cells<-rownames(annotations)
gene_expression<-gene_expression[, cells]


## 6. Heatmap with annotations
# 6.1 Set heatmap colors
library(RColorBrewer)
library("paletteer")

# 6.1.1 Cell types colors
celltype=unique(annotations$Cell_Subtype)
cb_palette<-paletteer_d("ggthemes::colorblind", 8)
cell_colors<-c(cb_palette, "#E41A1C", "#654CFFFF","#FFFFBFFF" )
names(cell_colors)<-celltype

# 6.1.2 Gene function colors
fun_colors<-brewer.pal(n = 4,name = "Paired")
names(fun_colors)<-c("Chaperone protein binding functions","Damage disease associated",
                     "Immune binding functions","Microglial state markers")

Function_mg <- c(Chaperone_pbf,Immune_bf,Microglial_sm,Damage_da)

# 6.1.2 Gene expression colors
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("deepskyblue2", "white", "red"))
col_fun(seq(-3, 3))


# 6.2 Set annotations
# 6.2.1 Column annotation
column_ha = HeatmapAnnotation(Cell_type =annotations$Cell_Subtype,
                                   Status = annotations$disease_ontology_label,
                                   col =  list(Status=c("control" = "cyan", "Parkinson's disease" = "salmon"),
                                               Cell_type = (df=c(cell_colors))),
                                   annotation_label=c("Cell type", "Status"),
                                   annotation_name_gp=gpar(fontsize = 10),
                                   annotation_legend_param = list(
                                     Status = list(
                                       title = "Status",
                                       labels = c("control", "Parkinson's disease"),
                                       nrow = 2,
                                       title_gp = gpar(fontsize = 9, fontface = "bold"),
                                       labels_gp = gpar(fontsize = 9)
                                     ),
                                     Cell_type = list(
                                       title = "Cell type",
                                       labels = names (cell_colors),
                                       annotation_label="Cell type",
                                       ncol = 3,
                                       title_gp = gpar(fontsize = 9, fontface = "bold"),
                                       labels_gp = gpar(fontsize = 7)
                                     )
                                   ))

# 6.2.2 Row annotation
row_ha = rowAnnotation(Molecular_function=Function_mg,
                          col=list(Molecular_function=c(fun_colors)),
                          show_annotation_name = FALSE,
                          annotation_legend_param = list(
                            Molecular_function = list(
                              title = "Molecular function",
                              labels = names (fun_colors),
                              title_gp = gpar(fontsize = 9, fontface = "bold"),
                              legend_gp = gpar(fun_colors),
                              labels_gp = gpar(fontsize = 9)
                            )))


# 6.3 Create heatmap
ht1_mg = Heatmap(gene_expression, name = "Scaled expression", col = col_fun,
                 row_names_side = "left",
                 cluster_rows = FALSE,
                 show_column_dend = FALSE, 
                 top_annotation = column_ha,
                 left_annotation =  row_ha,
                 column_order = order(annotations$disease_ontology_label,
                                      annotations$Cell_Type),
                 show_column_names = FALSE, use_raster= FALSE, 
                 row_title = NULL,
                 row_names_gp = gpar(fontsize = 7),
                 heatmap_legend_param=list(title_gp = gpar(fontsize = 9,fontface = "bold"),
                                           labels_gp=gpar(fontsize = 8))
)

# 6.4 Draw
draw(ht1_mg, heatmap_legend_side = "right", annotation_legend_side = "top", 
     align_annotation_legend="heatmap_center")

