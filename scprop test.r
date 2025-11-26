#prep library
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(clustree)
library("scProportionTest")

#read in objects
kamath_mg <- readRDS("kamath/kamath_mg.Rds")
kamath_astro <- readRDS("kamath/kamath_astro.Rds")

#proportionality test
#set the default to RNA!
DefaultAssay(kamath_mg) <- 'RNA'

# add the test
prop_test <- sc_utils(kamath_mg)

# Once the object is created, the permutation testing and bootstrapping can be run.
prop_test <- permutation_test(
  prop_test, cluster_identity = "Cell_Subtype",
  sample_1 = "normal", sample_2 = "Parkinson disease",
  sample_identity = "disease__ontology_label"
)
# Point-range plots of the results
#setting colors bc default is too light
pdf(paste("Kamath_mg proptest PD vs ctrl.pdf"), height = 5, width = 7)
permutation_plot(prop_test)+
  labs(title = "Kamath microglia subcluster proportions IPD vs control")+
	scale_color_manual(values = c("indianred3", "grey35"))
dev.off()
#repeat for kamath_astro

#going back to complete dataset for all cells proportionality
kamath_all <- readRDS("kamath/kamath_ann.Rds")
#remove unannotated cells
kamath_all <- subset(kamath_all, subset = Cell_Type == NaN, invert = TRUE)
#set the default to RNA!
DefaultAssay(kamath_all) <- 'RNA'

# add the test
prop_test <- sc_utils(kamath_all)

# Once the object is created, the permutation testing and bootstrapping can be run.
prop_test <- permutation_test(
  prop_test, cluster_identity = "Cell_Type",
  sample_1 = "normal", sample_2 = "Parkinson disease",
  sample_identity = "disease__ontology_label"
)

# Point-range plots of the results
pdf(paste("Kamath_all chi2 test proportionality PD vs ctrl.pdf"), height = 5, width = 7)
permutation_plot(prop_test)+
  labs(title = "Kamath cell type proportions IPD vs control")+
	scale_color_manual(values = c("indianred3", "grey35"))
dev.off()