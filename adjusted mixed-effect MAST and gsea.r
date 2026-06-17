
#template fr DGE analysis adjusted based on reviewer feedback
#library prep
library("ggplot2")
library("Seurat")
library("cowplot")
library("clustree")
library(patchwork)
library(dplyr)
library("Matrix")
library("MAST")
library(edgeR)

#read in objects
kamath_mg <- readRDS("kamath_mg.rds")

#check the following metadata slots exist
table(kamath_mg$disease__ontology_label)
length(unique(kamath_mg$donor_id))
table(kamath_mg$disease__ontology_label, kamath_mg$donor_id)

#expression profile for edger
mg_pb_counts <- AggregateExpression(
  kamath_mg,
  group.by = "donor_id",
  assays = "RNA",
  slot = "counts"
)$RNA
#strip leading "g"
colnames(mg_pb_counts) <- sub("^g", "", colnames(mg_pb_counts))
#build donor level metadata
pb_meta <- kamath_mg@meta.data[, c("donor_id", "disease__ontology_label", "sex")]
pb_meta <- pb_meta[!duplicated(pb_meta$donor_id), ]
rownames(pb_meta) <- pb_meta$donor_id
#align metadata to aggregated mtx
pb_meta <- pb_meta[colnames(mg_pb_counts), ]
#check the following is TRUE
identical(rownames(pb_meta), colnames(mg_pb_counts))
#comparing PD vs control
keep <- pb_meta$disease__ontology_label %in% c("Parkinson disease", "normal")
mg_pb_counts2 <- mg_pb_counts[, keep]
pb_meta2 <- pb_meta[keep, , drop = FALSE]

# set disease ontology as factor
pb_meta2$disease__ontology_label <- factor(pb_meta2$disease__ontology_label, levels = c("normal", "Parkinson disease"))
# Make sex a factor
pb_meta2$sex <- factor(pb_meta2$sex)

#edgeR wants integer coutns
mg_pb_counts2 <- round(mg_pb_counts2)
dge <- DGEList(counts = mg_pb_counts2)
#filter low-exp genes
keep_genes <- filterByExpr(dge, design = model.matrix(~ disease__ontology_label + sex, data = pb_meta2))
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
#TMM normalization
dge <- calcNormFactors(dge)
#design mtx
design <- model.matrix(~ disease__ontology_label + sex, data = pb_meta2)
#dispersion and QL fit
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
# PD vs normal test (coefficient name depends on factor coding)
colnames(design)
res <- glmQLFTest(fit, coef = "disease__ontology_labelParkinson disease")

#save results
mg_edger_tbl <- topTags(res, n = Inf)$table
mg_edger_tbl$gene <- rownames(mg_edger_tbl)
#add -log10 for later volcano plot
mg_edger_tbl$neglog10FDR <- -log10(mg_edger_tbl$FDR)
write.csv(mg_edger_tbl, "kamath_mg_pseudobulk_edgeR_PD_vs_normal.csv", row.names = FALSE)

#run enhanced volcano
dge_filtered <- mg_edger_tbl[mg_edger_tbl$PValue < 0.05, ]
nrow(dge_filtered)
#filter out the RP and CTD genes
dge_filtered_clean <- dge_filtered[
  !grepl("^(RPL|RPS|RP11-|CTD-)", dge_filtered$gene),
]

#filter top genes to look for labels of interest
top50_up_clean <- dge_filtered_clean[dge_filtered_clean$logFC > 0, ]
top50_up_clean <- top50_up_clean[order(top50_up_clean$logFC, decreasing = TRUE), ]
top50_up_clean <- head(top50_up_clean, 50)
top50_down_clean <- dge_filtered_clean[dge_filtered_clean$logFC < 0, ]
top50_down_clean <- top50_down_clean[order(top50_down_clean$logFC, decreasing = FALSE), ]
top50_down_clean <- head(top50_down_clean, 50)
write.csv(top50_up_clean,   "kamath_mg_top50_up_clean_logFC.csv", row.names = FALSE)
write.csv(top50_down_clean, "kamath_mg_top50_down_clean_logFC.csv", row.names = FALSE)
#list of labels based on above
volcano_labs <- c("CXCL13", "CCL18", "SLAMF7", "SLAMF8", "CLECL1", "CD52", "DCSTAMP", "LPL",
                  "CHORDC1", "NUPR1", "CES1", "SPR", "P4HA2", "TSKU", "AQP9", "COLEC11",
                  "CDKN2B", "GPNMB", "HSP90AB1", "LRRK2", "APOE", "HLA-DQA2", "HLA-DOB",
                  "IL2RA", "CH25H", "CD28", "CD244", "TNFRSF13C", "ITGB7", "OSM",
                  "CCL2", "SLC6A3", "SLC18A2", "EN1", "CALCR")

png("kamath_mg_edgeR_volcano_rawP.png",width = 8, height = 7, units = "in", res = 600)

EnhancedVolcano(
  mg_edger_tbl,
  lab = mg_edger_tbl$gene,
  selectLab = volcano_labs,
  x = 'logFC',
  y = 'PValue',
  title = 'Microglia pseudobulk (edgeR QL)',
  xlab = expression(Log[2]~fold~change),
  ylab = expression(-Log[10]~P~value),

  pCutoff = 0.1,
  FCcutoff = 0.5,

  pointSize = 2.5,
  labSize = 4,
  colAlpha = 0.8,

  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = Inf,

  ylim = c(0, 5)
)

dev.off()

#later i will try this with subtypes of mg but for now, let's try running GO analysis on this
library(clusterProfiler)
library(org.Hs.eg.db)
#read in saved results
mg_edger_tbl <- read.csv(
  "kamath_mg_pseudobulk_edgeR_PD_vs_normal.csv",
  stringsAsFactors = FALSE
)
#prep GO pipeline
geneList <- mg_edger_tbl$logFC
names(geneList) <- mg_edger_tbl$gene
geneList <- sort(geneList, decreasing = TRUE)

#run GO_BP
gsea_bp <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  eps           = 0
)

write.csv(as.data.frame(gsea_bp),
          "kamath_mg_GO_BP_GSEA_SYMBOL.csv",
          row.names = FALSE)
#run GO+MF
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

write.csv(as.data.frame(gsea_mf),
          "kamath_mg_GO_MF_GSEA_SYMBOL.csv",
          row.names = FALSE)

#visualize the results
library(enrichplot)
png("kamath_mg_gseaBP_dotplot.png",  width = 7, height = 5, units = "in", res = 600)
dotplot(gsea_bp, showCategory=10, font.size=7, split=".sign",
        label_format=40) + facet_grid(.~.sign) + 
  ggtitle("kamath MG gseaBP (edgeR pseudobulk contrast)") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold"),
        legend.text = element_text(size=7), legend.title = element_text(size=8))
dev.off()
#for MF
png("kamath_mg_gseaMF_dotplot.png", width = 7, height = 5, units = "in", res = 600)
dotplot(gsea_mf, showCategory=10, font.size=7, split=".sign",
        label_format=40) + facet_grid(.~.sign) + 
  ggtitle("kamath MG gseaMF (edgeR pseudobulk contrast)") +
  theme(plot.title = element_text(hjust = 0.4, face = "bold"),
        legend.text = element_text(size=7), legend.title = element_text(size=8))
dev.off()

#prepare for cnetplot adjusted for version ctrl issues
library(dplyr)
library(enrichplot)
library(ggplot2)
library(grid)   # unit()

# fold-change vector for coloring gene nodes (logFC)
fc <- mg_edger_tbl$logFC
names(fc) <- mg_edger_tbl$gene
fc <- fc[!is.na(names(fc)) & names(fc) != ""]
fc <- fc[!duplicated(names(fc))]

#pick top 5 activated and suppressed MF terms
mf_df <- as.data.frame(gsea_mf)

GO_up <- mf_df %>%
  filter(NES > 0) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  pull(Description)

GO_down <- mf_df %>%
  filter(NES < 0) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  pull(Description)

GO_cnet_terms <- c(GO_up, GO_down)

#compute termsim and plot cnet
x2 <- pairwise_termsim(gsea_mf)
#cnetplot layer with only category labels


png("kamath_mg_GSEA_MF_cnet.png", width = 10, height = 8, units = "in", res = 600)
print(p)
dev.off()

#go back and check MAST if sex is controlled for
kamath_mg <- readRDS("kamath_mg.rds")
Idents(kamath_mg) <- "disease__ontology_label"
DefaultAssay(kamath_mg) <- "RNA"
kamath_mg <- NormalizeData(kamath_mg)
#control for sex and donor_id
PDmarkers <- FindMarkers(kamath_mg, ident.1 = "Parkinson disease", ident.2 = "normal", test.use = "MAST", latent.vars = c("sex", "donor_id"))
PDmarkers_og <- FindMarkers(kamath_mg, ident.1 = "Parkinson disease", ident.2 = "normal", test.use = "MAST")
#prelim comparision of results
# Compare MAST results: with vs without covariates
# PDmarkers      = MAST + latent.vars (sex, donor_id)
# PDmarkers_og   = vanilla MAST

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# std and merge as dfs
prep_markers <- function(df, label) {
  df %>%
    rownames_to_column("gene") %>%
    mutate(model = label)
}

m_cov <- prep_markers(PDmarkers_sd, "MAST_cov_sex_donor")
m_og  <- prep_markers(PDmarkers_og, "MAST_vanilla")

#  detect p-value columns robustly across Seurat versions
# Common columns: p_val, p_val_adj, avg_log2FC / avg_logFC, pct.1, pct.2
required_cols <- c("gene")
stopifnot(all(required_cols %in% colnames(m_cov)), all(required_cols %in% colnames(m_og)))

# Identify the logFC column name
fc_col <- dplyr::case_when(
  "avg_log2FC" %in% colnames(m_cov) ~ "avg_log2FC",
  "avg_logFC"  %in% colnames(m_cov) ~ "avg_logFC",
  TRUE ~ NA_character_
)
if (is.na(fc_col)) stop("Could not find avg_log2FC/avg_logFC in PDmarkers output.")

p_col <- if ("p_val" %in% colnames(m_cov)) "p_val" else stop("Could not find p_val column.")
padj_col <- if ("p_val_adj" %in% colnames(m_cov)) "p_val_adj" else stop("Could not find p_val_adj column.")

# Merge by gene for side-by-side comparison
cmp <- m_cov %>%
  select(gene,
         !!fc_col := all_of(fc_col),
         p_val = all_of(p_col),
         p_val_adj = all_of(padj_col),
         pct.1 = any_of("pct.1"),
         pct.2 = any_of("pct.2")) %>%
  rename_with(~ paste0(., "_cov"), -gene) %>%
  left_join(
    m_og %>%
      select(gene,
             !!fc_col := all_of(fc_col),
             p_val = all_of(p_col),
             p_val_adj = all_of(padj_col),
             pct.1 = any_of("pct.1"),
             pct.2 = any_of("pct.2")) %>%
      rename_with(~ paste0(., "_og"), -gene),
    by = "gene"
  ) %>%
  mutate(
    neglog10_p_og  = -log10(p_val_og + 1e-300),
    neglog10_p_cov = -log10(p_val_cov + 1e-300),
    neglog10_fdr_og  = -log10(p_val_adj_og + 1e-300),
    neglog10_fdr_cov = -log10(p_val_adj_cov + 1e-300),
    d_logFC = .data[[paste0(fc_col, "_cov")]] - .data[[paste0(fc_col, "_og")]]
  )

# summary stats at a glance
sig_og_p   <- sum(cmp$p_val_og < 0.05, na.rm = TRUE)
sig_cov_p  <- sum(cmp$p_val_cov < 0.05, na.rm = TRUE)
sig_og_fdr  <- sum(cmp$p_val_adj_og < 0.05, na.rm = TRUE)
sig_cov_fdr <- sum(cmp$p_val_adj_cov < 0.05, na.rm = TRUE)

cat("\n=== Counts of significant genes ===\n")
cat("Vanilla MAST: p<0.05  =", sig_og_p,  "\n")
cat("MAST + cov:   p<0.05  =", sig_cov_p, "\n")
cat("Vanilla MAST: FDR<0.05=", sig_og_fdr, "\n")
cat("MAST + cov:   FDR<0.05=", sig_cov_fdr,"\n")

cat("\n=== Correlations (genes present in both) ===\n")
cat("logFC correlation:",
    cor(cmp[[paste0(fc_col, "_og")]], cmp[[paste0(fc_col, "_cov")]], use="pairwise.complete.obs"),
    "\n")
cat("-log10(p) correlation:",
    cor(cmp$neglog10_p_og, cmp$neglog10_p_cov, use="pairwise.complete.obs"),
    "\n")

# look for genes that changed most in significance
big_p_shift <- cmp %>%
  mutate(delta_neglog10p = neglog10_p_cov - neglog10_p_og) %>%
  arrange(desc(abs(delta_neglog10p))) %>%
  select(gene,
         starts_with("avg_log"), starts_with("p_val"), starts_with("p_val_adj"),
         delta_neglog10p) %>%
  head(30)

cat("\n=== Top 30 genes with biggest shift in -log10(p) (cov - og) ===\n")
print(big_p_shift)

#seems consistently the same logFC at least
#are there genes that changed directions?
signflip_genes <- cmp %>%
  filter(!is.na(avg_log2FC_og), !is.na(avg_log2FC_cov), sign(avg_log2FC_og) != sign(avg_log2FC_cov)) %>%
  select(gene,avg_log2FC_og, avg_log2FC_cov, p_val_og, p_val_cov, p_val_adj_og, p_val_adj_cov) %>%
  arrange(desc(abs(avg_log2FC_og - avg_log2FC_cov)))
dim(signflip_genes)
#good its empty

# overlap of top hits
topN <- 200
top_og <- m_og %>% arrange(p_val_adj) %>% slice_head(n = topN) %>% pull(gene)
top_cov <- m_cov %>% arrange(p_val_adj) %>% slice_head(n = topN) %>% pull(gene)

overlap <- length(intersect(top_og, top_cov))
cat("\n=== Overlap of top", topN, "by FDR ===\n")
cat("Overlap =", overlap, " (", round(100*overlap/topN, 1), "% of each list )\n", sep="")
#~20% overlap, remaining signature. proportional to edgeR results tbh

#save the overlap
preserved_DGE <- cmp %>%
  filter(p_val_adj_og < 0.05 & p_val_adj_cov < 0.05) %>%
  arrange(p_val_adj_cov)
write.csv(preserved_DGE,"kamath_mg_preservedDGEpostcorrection.csv",row.names = FALSE)

# diagnistic plots
# logFC concordance
png("MAST_change_in_effect_size.png", width=2200, height=2000, res=600)
print(
  ggplot(cmp, aes(x = .data[[paste0(fc_col, "_og")]],
                  y = .data[[paste0(fc_col, "_cov")]])) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(
      title = "MAST logFC: vanilla vs covariates",
      x = paste0("Vanilla ", fc_col),
      y = paste0("Covariates ", fc_col)
    )
)
dev.off()

# p-value inflation diagnostic: -log10(p) scatter
png("MAST_change_in_pinflation.png", width=2200, height=2000, res=300)
print(
  ggplot(cmp, aes(x = neglog10_p_og, y = neglog10_p_cov)) +
    geom_point(alpha = 0.25, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(
      title = "MAST significance: vanilla vs covariates",
      x = expression(-log[10](p)~"(vanilla)"),
      y = expression(-log[10](p)~"(covariates)")
    )
)
dev.off()

# Distribution of adjusted p-values
padj_long <- bind_rows(
  m_og %>% transmute(model="Vanilla", padj = .data[[padj_col]]),
  m_cov %>% transmute(model="Covariates", padj = .data[[padj_col]])
) %>% filter(!is.na(padj), padj <= 1, padj >= 0)

png("MAST_compare_FDR_hist.png", width=2200, height=1600, res=300)
print(
  ggplot(padj_long, aes(x = padj)) +
    geom_histogram(bins = 60) +
    facet_wrap(~model, ncol=1) +
    theme_bw() +
    labs(title="Adjusted p-value distributions (MAST)", x="FDR (p_val_adj)", y="Genes")
)
dev.off()

#save comparison table
write.csv(cmp, "MAST_compare_vanilla_vs_covariates_merged.csv", row.names = FALSE)
#also save both the MAST result
# original MAST
PDmarkers_og_out <- PDmarkers_og %>%
  rownames_to_column("gene")
write.csv(PDmarkers_og_out, "kamath_mg_MAST_og.csv",row.names = FALSE)

# MAST with covariates
PDmarkers_adj_out <- PDmarkers_sd %>%
  rownames_to_column("gene")
write.csv(PDmarkers_adj_out,"kamath_mg_MAST_adjusted.csv",row.names = FALSE)

#volcano and gsea done locally
#read in mixed effects MAST results
library(EnhancedVolcano)

mg_mast_adj <- read.csv(
  "C:/Users/belle/OneDrive - McGill University/Documents/Kamath_2026/mg_all_mixedMAST/kamath_mg_MAST_adjusted.csv",
  stringsAsFactors = FALSE
)

png("kamath_mg_volcano (mixed-MAST).png",width = 8, height = 8,units = "in",res = 600)

volcano_labs <- c("CXCL13", "CCL18", "SLAMF7", "SLAMF8", "CLECL1", "CD52", "DCSTAMP", "LPL",
                  "CHORDC1", "NUPR1", "CES1", "SPR", "P4HA2", "TSKU", "AQP9", "COLEC11",
                  "CDKN2B", "GPNMB", "HSP90AB1", "LRRK2", "APOE", "HLA-DQA2", "HLA-DOB",
                  "IL2RA", "CH25H", "CD28", "CD244", "TNFRSF13C", "ITGB7", "OSM",
                  "CCL2", "SLC6A3", "SLC18A2", "EN1", "CALCR")

EnhancedVolcano(
  mg_mast_adj,
  lab = mg_mast_adj$gene,
  selectLab = volcano_labs,
  x = "avg_log2FC",
  y = "p_val_adj",
  title = "Kamath microglia PD volcano (mixed-MAST)",
  subtitle = NULL,
  xlab = expression(Log[2]~fold~change),
  ylab = expression(-Log[10]~adjusted~P~value),
  
  pCutoff = 0.05,
  FCcutoff = 0.5,
  
  pointSize = 2.5,
  labSize = 4,
  colAlpha = 0.8,
  
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = Inf,
  
  legendPosition = "top",
  legendLabSize = 12,
  legendIconSize = 4.0,
  
  ylim = c(0, 200),
  xlim = c(-8, 8)
)
dev.off()

#repeat gsea BP & MF and dotplots
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# prep GO pipeline from mixed MAST results
geneList <- mg_mast_adj$avg_log2FC
names(geneList) <- mg_mast_adj$gene
geneList <- geneList[!is.na(names(geneList)) & names(geneList) != ""]
geneList <- geneList[!duplicated(names(geneList))]
geneList <- sort(geneList, decreasing = TRUE)

# run GO_BP
gsea_bp <- gseGO(
  geneList      = geneList,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  minGSSize     = 10,
  maxGSSize     = 500,
  pAdjustMethod = "BH",
  eps           = 0
)

write.csv(
  as.data.frame(gsea_bp),
  "kamath_mg_GO_BP_GSEA_SYMBOL (mixed-MAST).csv",
  row.names = FALSE
)

# run GO_MF
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
  "kamath_mg_GO_MF_GSEA_SYMBOL (mixed-MAST).csv",
  row.names = FALSE
)

# visualize GO_BP
png("kamath_mg_gseaBP_dotplot (mixed-MAST).png",
  width = 7, height = 5, units = "in", res = 600
)

dotplot(
  gsea_bp,
  showCategory = 10,
  font.size = 7,
  split = ".sign",
  label_format = 40
) +
  facet_grid(. ~ .sign) +
  ggtitle("kamath MG gseaBP (mixed-MAST)") +
  theme(
    plot.title = element_text(hjust = 0.4, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )

dev.off()

# visualize GO_MF
png(
  "kamath_mg_gseaMF_dotplot (mixed-MAST).png",
  width = 7, height = 5, units = "in", res = 600
)

dotplot(
  gsea_mf,
  showCategory = 10,
  font.size = 7,
  split = ".sign",
  label_format = 40
) +
  facet_grid(. ~ .sign) +
  ggtitle("kamath MG gseaMF (mixed-MAST)") +
  theme(
    plot.title = element_text(hjust = 0.4, face = "bold"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )

dev.off()