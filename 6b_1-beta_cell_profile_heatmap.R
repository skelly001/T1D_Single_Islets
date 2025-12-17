# ==============================================================================
# Beta Cell Profile Heatmap Visualization
# ==============================================================================
# Script: 6b_1-beta_cell_profile_heatmap.R
# Description: Creates comprehensive heatmap visualization of Beta Cell Profile
#              proteins across all donors using ComplexHeatmap. Features are centered
#              (zero-mean), and heatmaps include column annotation for both IIRS
#              (Islet Immune Response Signature) and BCP (Beta Cell Profile) centroid
#              values. Generates six donor-specific panels combined into a single figure.
#
# Input: - output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds
#        - output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds
#        - output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx
# Output: - output/RD6-beta_cell_profile/RD6b_1-beta_cell_profile_heatmap.png
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)
library(ComplexHeatmap)


# Import ------------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load protein abundance data with IIRS
m <- readRDS("output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds")

# Load IIRS data
m_imc <- readRDS(
	"output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds"
)

# Load BCP protein features
bcp_proteins <- openxlsx::read.xlsx(
	"output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx"
)


# Process -----------------------------------------------------------------

ht_opt$COLUMN_ANNO_PADDING = unit(3, "mm")
ht_opt$merge_legends = T
ht_opt$message = FALSE

# Subset to BCP proteins
m_sig <- m[match(bcp_proteins$gene_name, fData(m)$gene_name), ]

# Center features (zero-mean per protein)
exprs(m_sig) <- t(scale(t(exprs(m_sig)), center = T, scale = F))

# Calculate BCP centroid per sample
m_sig$bcp_centroid <- apply(exprs(m_sig), 2, \(x) {
	median(x, na.rm = T)
})

# Center IIRS for diverging color scale
m_imc$im_sig_centroid <- scale(m_imc$im_sig_centroid, center = T, scale = F)


## Generate Donor Heatmaps ----
donors <- unique(pData(m_sig)$donor)
h_hor <- vector("list", 6)
names(h_hor) <- donors


expr_col <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


for (i in donors) {
	d1 <- m_sig[, m_sig$donor == i]
	featureNames(d1) <- str_extract(featureNames(d1), "^[^\\.]+")
	d_IIRS <- m_imc[, m_imc$donor == i]

	# Column annotations for IIRS and BCP
	show_names <- i == "6539"
	ca1 <- HeatmapAnnotation(
		IIRS = d_IIRS$im_sig_centroid,
		col = list(IIRS = expr_col),
		show_annotation_name = show_names,
		simple_anno_size = unit(3.75, "mm")
	)

	ca2 <- HeatmapAnnotation(
		BCP = d1$bcp_centroid,
		col = list(BCP = expr_col),
		show_annotation_name = show_names,
		simple_anno_size = unit(3.75, "mm")
	)

	# Show legend only for last donor
	show_legend <- FALSE
	legend_param <- NULL

	if (i == "6539") {
		show_legend <- TRUE
		legend_param <- list(title = "Abundance")
	}

	h_hor[[i]] <- Heatmap(
		exprs(d1),
		top_annotation = c(ca1, ca2),
		cluster_rows = FALSE,
		clustering_distance_columns = "euclidean",
		clustering_method_columns = "average",
		show_column_names = F,
		column_title = i,
		show_heatmap_legend = show_legend,
		heatmap_legend_param = legend_param,
		col = expr_col
	)
}

png(
	filename = "output/RD6-beta_cell_profile/RD6b_1-beta_cell_profile_heatmap.png",
	height = 7.5,
	width = 10,
	units = "in",
	res = 600,
	type = "cairo"
)
draw(
	h_hor[[1]] + h_hor[[2]] + h_hor[[3]] + h_hor[[4]] + h_hor[[5]] + h_hor[[6]],
	gap = unit(3, "mm")
)
dev.off()
