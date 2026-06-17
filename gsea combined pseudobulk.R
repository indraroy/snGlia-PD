#combined pseudobulk analyses
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(EnhancedVolcano)

#read in the pseudobulk contrasts
combo_mg <- read.csv("C:/Users/belle/OneDrive - McGill University/scRNA seq/Kamath_Wang_Smajic/Wang_Smajc_Kamath_mg_pseudobulk_edgeR_PD_vs_normal.csv")
combo_astro <- read.csv("C:/Users/belle/OneDrive - McGill University/scRNA seq/Kamath_Wang_Smajic/Wang_Smajc_Kamath_astro_pseudobulk_edgeR_PD_vs_normal.csv")
#check contrast length

combo_astro_filt <- combo_astro %>%
  filter(
    FDR < 0.05,
    abs(logFC) > 0.5,
    logCPM > 0
  )

# Total number of genes passing filters
n_genes <- nrow(combo_astro_filt)

# Optional: split counts by direction
n_up <- combo_astro_filt %>% filter(logFC > 0) %>% nrow()
n_down <- combo_astro_filt %>% filter(logFC < 0) %>% nrow()

n_genes
n_up
n_down

#ok so that's not a lot, let's look at all of them
# Upregulated genes (logFC > 0)
up_genes <- combo_astro_filt %>%
  filter(logFC > 0) %>%
  pull(gene)

# Downregulated genes (logFC < 0)
down_genes <- combo_astro_filt %>%
  filter(logFC < 0) %>%
  pull(gene)

# Inspect
length(up_genes)
length(down_genes)

#see as strings
up_genes_str <- paste(up_genes, collapse = ", ")
down_genes_str <- paste(down_genes, collapse = ", ")
up_genes_str
down_genes_str

mg_up <- c("P2RX7","LGALS3", "LPL","PLIN2","SCARB2", "ALOX15B","GPNMB",  "LYZ","SNCA",
           "TRAF6", "IL27RA", "CLEC1A", "MIR155HG", "BCL10","MS4A14", "CYP27A1","ARG2", "SQSTM1", 
           "BNIP3L",  "DNAJB1", "HSPH1","HSPB1",  "CHORDC1", "STIP1", "HSPE1-MOB4", "P4HA2",
           "BACE2", "CRYZ", "HSPA4L", "CLIC2", "LUCAT1")
mg_down <- c("IFI44L", "RSAD2","MX2","OAS3", "TMEM173","CH25H","IL21R","P2RY6","LST1", "CIRBP",
             "FOLR2","PLD4", "HLA-DOA","TNFRSF21", "SNAP25", "IGF1", "ALDH1A1", "UNC13C", "KCNJ6")

#secondary lists for mg volcano that's just cytokine/immune receptor activity
mg_up <- c("P2RX7","IL27RA","CLEC1A", "TRAF6","BCL10", "GPNMB","MIR155HG")
mg_down <- c("TMEM173","IL21R","P2RY6","HLA-DOA","LST1","FOLR2","PLD4","TNFRSF21","CH25H")

#astrocytes
astro_up <- c( "SLC7A11", "GCLM", "NQO1", "SOD2", "LUCAT1", "FOXJ1", "PHGDH", "PSAT1", "SHMT2", "ASNS",
               "ATF4", "HERPUD1", "EDEM1", "HSPA5", "BAG3", "CHI3L1", "CD44", "SPARC", "IGFBP7", "GFPT2",
               "NUPR1", "EGLN3","NFKB2", "RELA", "IRF1", "CEBPG", "EIF4EBP1", "IL6ST", "NFKBIZ", "ABCA1",
               "MBOAT7",  "P4HA1", "PLOD2", "COL6A1", "CTGF","HSPE1-MOB4")
astro_down <- c("HES1", "HES5", "HEYL", "NRARP","SLC4A4", "SLC30A10","ADAMTSL1",
                "FGFR2", "FGFR3", "GFRA1", "PTN","CIRBP", "DDX60", "ACSL6", "ELOVL2",
                "GRIA1", "RBFOX1", "KCNJ6")
  
#test volcano
volcano_labs <- c(mg_up, mg_down)
volcano_labs2 <- c(astro_up, astro_down)

p<- EnhancedVolcano(
  combo_mg,
  lab = combo_mg$gene,
  selectLab = volcano_labs,
  x = "logFC",
  y = "FDR",
  title = "Combined microglia PD volcano (pseudobulk)",
  subtitle = NULL,
  xlab = expression(Log[2]~fold~change),
  ylab = expression(-Log[10]~adjusted~P~value),
  
  pCutoff = 0.05,
  FCcutoff = 0.5,
  ylim = c(0, 4),
  pointSize = 2.5,
  labSize = 4,
  colAlpha = 0.8,
  
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = Inf,
  
  legendPosition = "top",
  legendLabSize = 12,
  legendIconSize = 4.0
)

#save volcano
ggsave(
  filename = "combined_microglia_volcano_immune.png",
  plot = p,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 600
)

#volcano for astrocytes
p <- EnhancedVolcano(
  combo_astro,
  lab = combo_astro$gene,
  selectLab = volcano_labs2,
  x = "logFC",
  y = "FDR",
  title = "Combined astrocyte PD volcano (pseudobulk)",
  subtitle = NULL,
  xlab = expression(Log[2]~fold~change),
  ylab = expression(-Log[10]~FDR),
  
  pCutoff = 0.05,
  FCcutoff = 0.5,
  ylim = c(0, 6),
  
  pointSize = 2.5,
  labSize = 4,
  colAlpha = 0.8,
  
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = Inf,
  
  legendPosition = "top",
  legendLabSize = 12,
  legendIconSize = 4.0
)

#save volcano
ggsave(
  filename = "combined_astro_volcano.png",
  plot = p,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 600
)

#GO_MF analysis of both contrasts
library(clusterProfiler)
library(org.Hs.eg.db)

#run for mg first, later will be replaced by astro
geneList <- combo_mg$logFC
names(geneList) <- combo_mg$gene
geneList <- sort(geneList, decreasing = TRUE)
gsea_mf <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "MF",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  eps           = 0
)
write.csv(
  as.data.frame(gsea_mf),
  "combined_pseudobulk_microglia_GO_MF_GSEA_SYMBOL.csv",
  row.names = FALSE
)

#simplified dotplot bc lotta redundancy
# Simplify redundant GO terms
gsea_mf_simplified <- simplify(
  gsea_mf,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
# Make simplified dotplot
p_simple <- dotplot(
  gsea_mf_simplified,
  showCategory = 10,
  font.size = 7,
  split = ".sign",
  label_format = 40,
  orderBy = "NES"
) +
  facet_grid(. ~ .sign) +
  ggtitle("GO_molecular function term enrichment, PD microglia (multi-dataset)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
#save plot
ggsave(
  filename = "combined_pseudobulk_microglia_GO_MF_dotplot_simplified.png",
  plot = p_simple,
  width = 7.5,
  height = 5.5,
  units = "in",
  dpi = 600
)

#repeat for astrocytes
geneList <- combo_astro$logFC
names(geneList) <- combo_astro$gene
geneList <- sort(geneList, decreasing = TRUE)
gsea_mf <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "MF",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  eps           = 0
)
write.csv(
  as.data.frame(gsea_mf),
  "combined_pseudobulk_astrocytes_GO_MF_GSEA_SYMBOL.csv",
  row.names = FALSE
)
gsea_mf_simplified <- simplify(
  gsea_mf,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)
p_simple <- dotplot(
  gsea_mf_simplified,
  showCategory = 10,
  font.size = 7,
  split = ".sign",
  label_format = 40,
  orderBy = "NES"
) +
  facet_grid(. ~ .sign) +
  ggtitle("GO molecular function term enrichment, PD astrocytes (multi-dataset)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
p_simple
ggsave(
  filename = "C:/Users/belle/OneDrive - McGill University/scRNA seq/Kamath_Wang_Smajic/combined_pseudobulk_astrocytes_GO_MF_dotplot_simplified2.png",
  plot = p_simple,
  width = 9.5,
  height = 7.5,
  units = "in",
  dpi = 600
)

#check overlap between the two df
astro_filt <- combo_astro %>%
  filter(
    FDR < 0.05,
    abs(logFC) > 0.5,
    logCPM > 0
  ) %>%
  dplyr::select(gene, logFC, FDR)

mg_filt <- combo_mg %>%
  filter(
    FDR < 0.05,
    abs(logFC) > 0.5,
    logCPM > 0
  ) %>%
  dplyr::select(gene, logFC, FDR)

#overlap
overlap <- inner_join(
  astro_filt,
  mg_filt,
  by = "gene",
  suffix = c("_astro", "_mg")
)
#directionality
overlap <- overlap %>%
  mutate(
    direction_astro = ifelse(logFC_astro > 0, "Up", "Down"),
    direction_mg    = ifelse(logFC_mg > 0, "Up", "Down"),
    same_direction  = direction_astro == direction_mg
  )
table(overlap$same_direction)
#alright these are all true
#looking for up vs down
overlap_up <- overlap %>%
  filter(logFC_astro > 0 & logFC_mg > 0)
overlap_down <- overlap %>%
  filter(logFC_astro < 0 & logFC_mg < 0)
#get genes as strings
overlap_up_genes <- overlap_up$gene
overlap_down_genes <- overlap_down$gene
overlap_up_str <- paste(overlap_up_genes, collapse = ", ")
overlap_down_str <- paste(overlap_down_genes, collapse = ", ")
overlap_up_str
overlap_down_str

#scatterplot of the overlapping genes' connection
library(dplyr)
library(ggplot2)
library(ggrepel)

# Quadrant classification (only concordant)
overlap <- overlap %>%
  mutate(
    quadrant = case_when(
      logFC_astro > 0 & logFC_mg > 0 ~ "Up in both",
      logFC_astro < 0 & logFC_mg < 0 ~ "Down in both"
    )
  )

cor_val <- cor(overlap$logFC_astro, overlap$logFC_mg)

p_scatter <- ggplot(overlap, aes(x = logFC_astro, y = logFC_mg)) +
  
  # Points
  geom_point(aes(color = quadrant), alpha = 0.7, size = 2.3) +
  
  # Reference lines
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  
  # Diagonal
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dotted") +
  
  # Labels (black)
  geom_text_repel(
    data = overlap %>% filter(gene %in% genes_to_label),
    aes(label = gene),
    color = "black",
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  
  # Correlation
  annotate(
    "text",
    x = min(overlap$logFC_astro),
    y = max(overlap$logFC_mg),
    label = paste0("r = ", round(cor_val, 2)),
    hjust = 0,
    vjust = 1,
    size = 5
  ) +
  
  # Colors
  scale_color_manual(
    values = c(
      "Up in both" = "#D55E00",
      "Down in both" = "#0072B2"
    )
  ) +
  
  labs(
    title = "Concordant gene expression, microglia & astrocytes",
    x = expression(Astrocyte~Log[2]~fold~change),
    y = expression(Microglia~Log[2]~fold~change)
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    plot.title = element_text(
      size = 11,
      hjust = 0,
      face = "bold"
    ),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
  ) +
  
  coord_equal()

p_scatter
#save for publication
ggsave(
  "astro_vs_microglia_scatter.png",
  plot = p_scatter,
  width = 6,
  height = 6,   # 👈 make equal
  units = "in",
  dpi = 600
)