# ==============================================================================
# Miscellaneous Quality Control Plots for Islet Markers
# ==============================================================================
# Script: 10a_1-misc_QC_plots.R
# Description: Generates quality control visualizations for key islet hormone
#              markers and protein completeness. Creates boxplots showing insulin
#              (INS) and glucagon (GCG) intensity distributions by donor, scatter
#              plot examining INS vs GCG correlation, and barplots summarizing
#              observed protein counts across donors with >50% completeness.
#
# Input: - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
# Output: - output/RD10-misc/RD10a_1-INS_intensity_boxplot.png
#         - output/RD10-misc/RD10a_1-GCG_intensity_boxplot.png
#         - output/RD10-misc/RD10a_1-INS_vs_GCG_scatterplot.png
#         - output/RD10-misc/RD10a_1-observed_proteins_QC_barplot.png
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)


# Import ------------------------------------------------------------------

# Load protein abundance data with centroids
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")


# INS & GCG QC Plots ------------------------------------------------------

dir.create("output/RD10-misc")


# Define marker colors
rcartocolor::carto_pal(11, "Safe")
ins_color <- "#CC6677"
gcg_color <- "#88CCEE"

# Helper function: protein intensity distribution by donor
box_plot <- \(protein_abundance, color) {
	ggplot(pData(m), aes(x = donor, y = {{ protein_abundance }})) +
		geom_boxplot(
			fill = color,
			outlier.shape = 21,
			outlier.alpha = 0.6,
			outlier.size = 1.5
		) +
		theme_bw() +
		theme(
			panel.grid.major.x = element_blank(),
			axis.ticks.x = element_blank(),
			axis.text.x = element_text(
				angle = 90,
				hjust = 1,
				vjust = 0.5,
				color = "black"
			),
			axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			legend.position = "none"
		)
}

## INS Intensity ----
INS_intensity <- exprs(m)["INS", ]
box_plot(INS_intensity, ins_color) +
	labs(x = NULL, y = "INS")

ggsave(
	filename = "output/RD10-misc/RD10a_1-INS_intensity_boxplot.png",
	height = 4,
	width = 4,
	dpi = 600
)


## GCG Intensity ----
GCG_intensity <- exprs(m)["GCG", ]
box_plot(GCG_intensity, gcg_color) +
	labs(x = NULL, y = "GCG")

ggsave(
	filename = "output/RD10-misc/RD10a_1-GCG_intensity_boxplot.png",
	height = 4,
	width = 4,
	dpi = 600
)


## INS vs GCG Correlation ----
ggplot(pData(m), aes(x = INS_intensity, y = GCG_intensity)) +
	geom_point(
		alpha = 0.5,
		size = 1.5
	) +
	labs(x = "INS", y = "GCG") +
	theme_bw() +
	theme(
		panel.grid.minor.x = element_blank(),
		axis.text.x = element_text(color = "black"),
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 14),
		legend.position = "none"
	) +
	coord_fixed(ratio = 1) +
	xlim(c(-4.4, 3)) +
	ylim(c(-4.4, 3))

ggsave(
	filename = "output/RD10-misc/RD10a_1-INS_vs_GCG_scatterplot.png",
	height = 4,
	width = 4,
	dpi = 600
)


# Protein Count per Donor -------------------------------------------------

res <- data.frame(
	Donor = unique(pData(m)$donor),
	Class = str_to_sentence(unique(pData(m)[, c("donor", "class")])$class),
	Donor = paste0(c(rep('PreT1D-', 3), rep('Ctrl-', 3)), rep(1:3, 2)),
	num_prot_all = NA,
	num_prot_50pct = NA,
	num_prot_75pct = NA
) %>%
	mutate(
		Class = if_else(Class == "Case", "mAAb+", "Non-diabetic"),
		Class = factor(Class, levels = c("mAAb+", "Non-diabetic"))
	)


for (i in unique(pData(m)$donor)) {
	m_lp <- m[, pData(m)$donor == i]

	res[res$Donor == i, "num_prot_all"] <- nrow(exprs(m_lp))

	m_lp <- filter_by_occurrence(m_lp, 0.5)
	res[res$Donor == i, "num_prot_50pct"] <- nrow(exprs(m_lp))

	m_lp <- filter_by_occurrence(m_lp, 0.75)
	res[res$Donor == i, "num_prot_75pct"] <- nrow(exprs(m_lp))
}


p1 <- ggplot(res, aes(fct_inorder(Donor), num_prot_50pct, fill = Class)) +
	geom_bar(stat = "identity") +
	ylab("Number of Proteins") +
	xlab("Donor") +
	theme_classic() +
	scale_y_continuous(breaks = c(seq(0, 6000, 1000))) +
	theme(
		axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 10),
		axis.text.y = element_text(size = 10),
		axis.title = element_text(size = 12),
		legend.text = element_text(size = 10),
		legend.title = element_text(size = 11),
		legend.key.size = unit(0.15, "in"),
		legend.position = "bottom"
	) +
	guides(fill = guide_legend(nrow = 2))


ggsave(
	filename = "output/RD10-misc/RD10a_1-observed_proteins_QC_barplot.png",
	plot = p1,
	width = 2.1,
	height = 3,
	dpi = 600
)
