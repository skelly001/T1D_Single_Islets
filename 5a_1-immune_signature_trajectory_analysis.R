# ==============================================================================
# IIRS Trajectory Analysis via UMAP Dimensionality Reduction
# ==============================================================================
# Script: 5a_1-immune_signature_trajectory_analysis.R
# Description: Performs UMAP dimensionality reduction on the Islet Immune Response
#              Signature (IIRS) protein data to visualize trajectory patterns across
#              single islets. Generates UMAP plots showing IIRS centroid values for
#              all donors combined and separately by individual donor to reveal
#              immune response progression patterns in mAAb+ and control donors.
#
# Input: - output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds
#        - output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds
# Output: - output/RD5-islet_immune_response_signature_analysis/RD5a_1-UMAP_IIRS_trajectory.png
#         - output/RD5-islet_immune_response_signature_analysis/RD5a_1-UMAP_IIRS_trajectory_donor_faceted.png
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

# Load IIRS data with case donor adjustment
m_imc_noimpute <- readRDS(
	"output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds"
)


# Process -----------------------------------------------------------------

# Impute missing values for UMAP
m_imc <- impute_msnset(m_imc_noimpute)


# Trajectory analysis -----------------------------------------------------

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
	mutate(im_sig_centroid = im_sig_cor, donor = pData(m_imc)$donor) %>%
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

# Combined UMAP (all donors)
im_sig_limit <- max(abs(range(rd2_plot$im_sig_centroid, na.rm = T)))
p_main <- ggplot(rd2_plot, aes(UMAP1, UMAP2, color = im_sig_centroid)) +
	geom_point(size = 3) +
	labs(color = "IIRS") +
	plot_color(im_sig_limit) +
	plot_theme

ggsave(
	"output/RD5-islet_immune_response_signature_analysis/RD5a_1-UMAP_IIRS_trajectory.png",
	plot = p_main,
	width = 9,
	height = 6,
	dpi = 600
)


# Faceted UMAP by donor
im_sig_limit <- max(abs(range(rd2_plot$im_sig_centroid, na.rm = T)))
p1 <- ggplot(rd2_plot, aes(UMAP1, UMAP2, color = im_sig_centroid)) +
	geom_point() +
	facet_wrap("donor") +
	labs(color = "IIRS") +
	plot_color(im_sig_limit) +
	plot_theme

ggsave(
	"output/RD5-islet_immune_response_signature_analysis/RD5a_1-UMAP_IIRS_trajectory_donor_faceted.png",
	plot = p1,
	width = 9,
	height = 4,
	dpi = 600
)

# Overlay tissue section
p_lobe <- pData(m_imc)[, c("sample_name", "Block")]

rd2_plot <- rd2_plot %>%
	left_join(p_lobe, by = c("rowname" = "sample_name")) %>%
	mutate(lobe = str_extract(Block, "[[:alpha:]]{2}$"))

p2 <- ggplot(rd2_plot, aes(UMAP1, UMAP2, color = lobe)) +
	geom_point() +
	facet_wrap("donor") +
	labs(color = "IIRS") +
	plot_theme
p2

# IIRS by tissue section boxplot
im_df <- as.data.frame.MSnSet(m_imc) %>%
	mutate(lobe = str_extract(Block, "[[:alpha:]]{2}$"))

ggplot(im_df, aes(lobe, im_sig_centroid)) +
	geom_boxplot() +
	facet_wrap("donor")
