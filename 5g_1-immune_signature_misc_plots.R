# ==============================================================================
# Immune Signature Miscellaneous Plots
# ==============================================================================
# Script: 5g_1-immune_signature_misc_plots.R
# Description: Generates supplementary visualizations showing IIRS distribution,
#              protein-protein correlations, and co-expression patterns. Creates
#              (1) intra-donor variance boxplots, (2) STAT1 vs WARS1 correlation
#              scatterplots comparing case and control donors, and (3) IIRS protein
#              correlation heatmaps using Pearson correlation coefficients.
#
# Input: - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
#        - output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds
#        - output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds
# Output: - output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_intradonor_variance_boxplot.png
#         - output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_STAT1_WARS1_correlation_scatterplots.png
#         - output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_correlation_heatmap.png
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(corrr)


# Import -----------------------------------------------------------------

# Load protein abundance data with centroids
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")

# Load IIRS data
m_imc <- readRDS(
	"output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds"
)

# Load IIRS protein features
im_sig <- readRDS("output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds")


# IIRS Intra-Donor Variance ---------------------------------------------------

var_plot_data <- pData(m_imc)[, c("class", "donor", "im_sig_centroid")] %>%
	mutate(
		Class = if_else(class == "case", "mAAb+", "Non-Diabetic"),
		Class = factor(Class, levels = c("mAAb+", "Non-Diabetic"))
	)

var_plot <- ggplot(var_plot_data, aes(donor, im_sig_centroid, fill = Class)) +
	geom_boxplot(outliers = T, width = 0.4, linewidth = 0.2) +
	theme_bw() +
	xlab("Donor") +
	ylab("Islet Immune Response Signature (IIRS)") +
	scale_fill_manual(values = c("#ff0000", "#2222ff")) +
	theme(
		axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
		legend.position = "bottom"
	) +
	guides(fill = guide_legend(nrow = 2))

ggsave(
	"output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_intradonor_variance_boxplot.png",
	plot = var_plot,
	width = 1.9,
	height = 6.5,
	units = "in",
	dpi = 600
)


# IIRS Protein Correlation Heatmap --------------------------------------------

# Subset to IIRS proteins
m_im <- m[fData(m)$UniProtAcc %in% im_sig$UniProtAcc, ]

donor_case <- 6521
donor_cont <- 6178

# Helper function: generate correlation heatmap for a donor
make_cor_heatmap <- function(donor_id, title, show_legend = FALSE) {
	correlate(t(exprs(m_im)[sort(featureNames(m_im)), m_im$donor == donor_id])) %>%
		shave() %>%
		pivot_longer(cols = !matches("term"), names_to = "protein2", values_to = "cor") %>%
		dplyr::rename(protein1 = term) %>%
		drop_na(cor) %>%
		ggplot(aes(protein1, protein2, fill = cor)) +
		geom_tile() +
		scale_fill_gradientn(colors = c("blue", "#FFFFFF", "red"), limits = c(-1, 1)) +
		theme_classic() +
		ggtitle(title) +
		theme(
			axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
			axis.text.y = element_text(size = 6),
			axis.title = element_blank(),
			plot.title = element_text(hjust = 0.5, size = 12),
			legend.position = if (show_legend) "right" else "none",
			legend.title = element_text(size = 10),
			legend.text = element_text(size = 8),
			legend.key.size = unit(0.3, "cm")
		) +
		labs(fill = "Pearson\nCorrelation")
}

c1 <- make_cor_heatmap(donor_case, str_glue("mAAb+\n{donor_case}"))
c2 <- make_cor_heatmap(donor_cont, str_glue("Non-Diabetic\n{donor_cont}"), show_legend = TRUE)

c3 <- c1 + c2

ggsave(
	"output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_correlation_heatmap.png",
	plot = c3,
	width = 8.5,
	height = 4.2,
	dpi = 600
)
