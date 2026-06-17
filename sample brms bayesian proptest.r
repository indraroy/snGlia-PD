#brms as an alternative bayesian metho for scproptest
#compositional analysis packages
library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(readr)
library(brms)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(forcats)
library(scales)
#sanity check
cmdstanr::cmdstan_version()

#read in object
kamath_all <- readRDS("kamath_all.rds")
#check metadata
colnames(kamath_all@meta.data)

#pull down metadata for PD vsnormal
meta <- kamath_all@meta.data %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell_barcode")
required_cols <- c("donor_id", "Cell_Type", "disease__ontology_label")
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}
meta2 <- meta %>%
  filter(disease__ontology_label %in% c("Parkinson disease", "normal")) %>%
  mutate(
    condition = case_when(
      disease__ontology_label == "normal" ~ "CTRL",
      disease__ontology_label == "Parkinson disease" ~ "PD",
      TRUE ~ NA_character_
    ),
    condition = factor(condition, levels = c("CTRL", "PD")),
    donor_id = as.character(donor_id),
    Cell_Type = as.character(Cell_Type)
  ) %>%
  filter(!is.na(condition), !is.na(donor_id), !is.na(Cell_Type), Cell_Type != "")

cat("Cells retained:", nrow(meta2), "\n")
cat("Unique donors retained:", dplyr::n_distinct(meta2$donor_id), "\n")
table(meta2$condition)

#build donor-level cell count table
donor_celltype_counts_long <- meta2 %>%
  count(donor_id, condition, Cell_Type, name = "n_cells")

donor_celltype_counts_wide <- donor_celltype_counts_long %>%
  tidyr::pivot_wider(
    names_from = Cell_Type,
    values_from = n_cells,
    values_fill = 0
  ) %>%
  ungroup()

#identify celltype col
celltype_cols <- setdiff(colnames(donor_celltype_counts_wide), c("donor_id", "condition"))

# Drop donors with zero total cells across retained cell types, just in case
donor_celltype_counts_wide <- donor_celltype_counts_wide %>%
  mutate(total_cells = rowSums(across(all_of(celltype_cols)))) %>%
  filter(total_cells > 0)

# Quick check
print(donor_celltype_counts_wide[, c("donor_id", "condition", "total_cells")])
cat("Cell type columns:\n")
print(celltype_cols)
#set olig as reference
refcat <- "olig"
#build multinomial response matrix
model_df <- donor_celltype_counts_wide

# Put reference category first for readability / consistency
ordered_celltype_cols <- c(refcat, setdiff(celltype_cols, refcat))

# brms expects a matrix response for multinomial counts
count_mat <- as.matrix(model_df[, ordered_celltype_cols, drop = FALSE])

# Safety checks
if (!is.numeric(count_mat)) {
  storage.mode(count_mat) <- "numeric"
}
if (ncol(count_mat) < 2) {
  stop("Need at least 2 cell-type categories in the count matrix.")
}
if (any(rowSums(count_mat) <= 0)) {
  stop("At least one donor has zero total cells after filtering.")
}

model_df$total <- rowSums(count_mat)
model_df$counts <- count_mat

str(model_df$counts)
dim(model_df$counts)
colnames(model_df$counts)

#fit brms model
fit_kamath_all <- brm(
  formula = counts | trials(total) ~ condition,
  data = model_df,
  family = multinomial(refcat = refcat),
  chains = 4,
  cores = 4,
  iter = 4000,
  warmup = 2000,
  seed = 1234,
  backend = "cmdstanr",
  file = "kamath_all_celltype_composition_brms"
)

saveRDS(model_df, "kamath_all_celltype_composition_input_df.rds")
saveRDS(fit_kamath_all, "kamath_all_celltype_composition_brmsfit.rds")

#get summaries
print(summary(fit_kamath_all))

#getting posterior predictions from my model
newdata <- data.frame(
  condition = factor(c("CTRL", "PD"), levels = c("CTRL", "PD")),
  total = 1
)

epred <- posterior_epred(
  fit_kamath_all,
  newdata = newdata,
  re_formula = NA
)
n_draws <- dim(epred)[1]
cell_types <- colnames(model_df$counts)
#convert to tiduform
epred_df <- purrr::map_dfr(1:n_draws, function(d) {
  tibble(
    draw = d,
    condition = rep(c("CTRL", "PD"), each = length(cell_types)),
    Cell_Type = rep(cell_types, times = 2),
    prob = c(epred[d, 1, ], epred[d, 2, ])
  )
})

#compute log2fc per posterior draw
fc_draws <- epred_df %>%
  pivot_wider(names_from = condition, values_from = prob) %>%
  mutate(
    log2_fc = log2(PD / CTRL)
  )
#summarize posterior distribution
fc_summary <- fc_draws %>%
  group_by(Cell_Type) %>%
  summarise(
    mean_log2_fc = mean(log2_fc),
    lower = quantile(log2_fc, 0.025),
    upper = quantile(log2_fc, 0.975),
    prob_positive = mean(log2_fc > 0),
    .groups = "drop"
  ) %>%
  arrange(mean_log2_fc)
#reorder
fc_summary$Cell_Type <- factor(
  fc_summary$Cell_Type,
  levels = rev(unique(fc_summary$Cell_Type))
)

#plot (as dot and whiskers now)
p_fc <- ggplot(fc_summary,
               aes(y = Cell_Type, x = mean_log2_fc)) +
  
  # reference lines
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = c(-0.56, 0.56), linetype = "dotted") +
  
  # point estimate
  geom_point(size = 3, color = "#1f3b73") +
  
  # 95% CI whiskers
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2,
    color = "#1f3b73"
  ) +
  
  # annotation
  annotate(
    "text",
    x = min(fc_summary$lower) * 0.9,
    y = Inf,
    label = "Whiskers = 95% credible interval (posterior)",
    hjust = 0,
    vjust = 1.5,
    size = 4
  ) +
  
  coord_cartesian(clip = "off") +
  
  theme_bw(base_size = 14) +
  theme(
    plot.margin = margin(10, 20, 10, 10)
  ) +
  
  labs(
    title = "Log2FC (PD vs control) Bayesian model",
    x = "Model log2FC (posterior, Parkinson's disease vs control)",
    y = "Cell Type"
  )

p_fc
#export high res
ggsave(filename = "Cell proportionality graphs/kamath_all_brms_log2FC_plot.png",
       plot = p_fc,
       width = 8,
       height = 6,
       units = "in",
       dpi = 600
)