# ==============================================================================
# Beta Cell Profile Intra-Donor Variance Visualization
# ==============================================================================
# Script: 7e_1-beta_cell_profile_misc_plots.R
# Description: Generates supplementary visualization showing Beta Cell Profile
#              (BCP) centroid distribution across donors and disease status. Creates
#              boxplots stratified by donor and disease class (mAAb+ vs
#              Non-Diabetic) to illustrate intra-donor variance in beta cell function.
#
# Input: - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
#        - output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx
# Output: - output/RD7-beta_cell_profile_analysis/RD7e_1-BCP_intradonor_variance_boxplot.png
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)


# Import -----------------------------------------------------------------

# Load protein abundance data with centroids
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")

# Load BCP protein features
bcp_proteins <- openxlsx::read.xlsx(
	"output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx"
)


# BCP Intra-Donor Variance ------------------------------------------------

m_bcp <- m[fData(m)$gene_name %in% bcp_proteins$gene_name, ]

var_plot_data <- pData(m_bcp)[, c("class", "donor", "bcp_centroid")] %>%
	mutate(
		Class = if_else(class == "case", "mAAb+", "Non-Diabetic"),
		Class = factor(Class, levels = c("mAAb+", "Non-Diabetic"))
	)

var_plot <- ggplot(var_plot_data, aes(donor, bcp_centroid, fill = Class)) +
	geom_boxplot(outliers = T, width = 0.4, linewidth = 0.2) +
	theme_bw() +
	xlab("Donor") +
	ylab("Beta Cell Profile (BCP)") +
	scale_fill_manual(values = c("#ff0000", "#2222ff")) +
	theme(
		axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
		legend.position = "bottom"
	) +
	guides(fill = guide_legend(nrow = 2))

ggsave(
	"output/RD7-beta_cell_profile_analysis/RD7e_1-BCP_intradonor_variance_boxplot.png",
	plot = var_plot,
	width = 1.9,
	height = 6.5,
	dpi = 600
)
