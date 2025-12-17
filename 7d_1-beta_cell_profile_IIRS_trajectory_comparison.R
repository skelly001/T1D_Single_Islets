# ==============================================================================
# Beta Cell Profile and IIRS Trajectory Comparison Using UMAP
# ==============================================================================
# Script: 7d_1-beta_cell_profile_IIRS_trajectory_comparison.R
# Description: Compares Islet Immune Response Signature (IIRS) and Beta Cell
#              Profile (BCP) trajectory patterns using UMAP dimensionality reduction.
#              Performs imputation, calculates BCP centroid, and generates faceted
#              UMAP plots overlaying both IIRS and BCP centroid values in the IIRS UMAP
#              space to visualize the relationship between immune response and
#              beta cell function across donors.
#
# Input: - output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds
#        - output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds
#        - output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx
# Output: - output/RD7-beta_cell_profile_analysis/RD7d_1-UMAP_IIRS_BCP_trajectory_comparison_donor_faceted.png
# ==============================================================================

library(SingleCellExperiment)
library(MSnSet.utils)
library(patchwork)
library(uwot)
library(tidyverse)


# Import ------------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load batch-corrected protein abundance data
m <- readRDS("output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds")

# Load IIRS data
m_imc_noimpute <- readRDS(
	"output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds"
)

# Load BCP protein features
bcp_proteins <- openxlsx::read.xlsx(
	"output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx"
)


# Process -----------------------------------------------------------------

# Impute missing values for UMAP
m_imc <- impute_msnset(m_imc_noimpute)

# Subset to BCP proteins
m_bcp <- m[fData(m)$gene_name %in% bcp_proteins$gene_name, ]

# Calculate BCP centroid per sample
pData(m_bcp)$bcp_centroid <- apply(exprs(m_bcp), 2, \(x) {
	median(x, na.rm = TRUE)
})


# Trajectory comparison -----------------------------------------------------

## Dimensionality Reduction ----

# Compute UMAP embedding
rd2 <- umap(t(exprs(m_imc)), scale = T, seed = 0)
colnames(rd2) <- c('UMAP1', 'UMAP2')
identical(rownames(rd2), sampleNames(m_imc))


# Center IIRS values for diverging color scale
im_sig <- pData(m_imc)$im_sig_centroid
im_sig_range <- range(im_sig, na.rm = T)
cor_factor <- max(im_sig, na.rm = T) - (im_sig_range[[2]] - im_sig_range[[1]]) / 2
im_sig_cor <- im_sig - cor_factor

# Prepare data for visualization
rd2_plot <- rd2 %>%
	as.data.frame() %>%
	mutate(
		im_sig_centroid = im_sig_cor,
		bcp_centroid = pData(m_bcp)$bcp_centroid,
		donor = pData(m_bcp)$donor
	) %>%
	rownames_to_column()


## Visualization ----

### Plot Components ----

# Color scale function
plot_color <- function(limit) {
	paletteer::scale_color_paletteer_c(
		"grDevices::Zissou 1",
		direction = 1,
		limits = c(-limit, limit)
	)
}

plot_theme <- theme_bw()


### Generate Plots ----

# Faceted UMAP by donor
im_sig_limit <- max(abs(range(rd2_plot$im_sig_centroid, na.rm = T)))
p1 <- ggplot(rd2_plot, aes(UMAP1, UMAP2, color = im_sig_centroid)) +
	geom_point() +
	facet_wrap("donor") +
	labs(color = "IIRS") +
	plot_color(im_sig_limit) +
	plot_theme

bcp_limit <- max(abs(range(rd2_plot$bcp_centroid, na.rm = T)))
p2 <- ggplot(rd2_plot, aes(UMAP1, UMAP2, color = bcp_centroid)) +
	geom_point() +
	facet_wrap("donor") +
	labs(color = "BCP") +
	plot_color(bcp_limit) +
	plot_theme

(p3 <- p1 / p2)

ggsave(
	"output/RD7-beta_cell_profile_analysis/RD7d_1-UMAP_IIRS_BCP_trajectory_comparison_donor_faceted.png",
	plot = p3,
	width = 9,
	height = 8,
	dpi = 600
)
