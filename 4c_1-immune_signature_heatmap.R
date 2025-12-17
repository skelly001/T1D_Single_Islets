# ==============================================================================
# Immune Signature Heatmap Visualization
# ==============================================================================
# Script: 4c_1-immune_signature_heatmap.R
# Description: Generates heatmap visualizations of the Islet Immune Response
#              Signature (IIRS) proteins across all samples. Creates separate
#              heatmaps for different sample classes (INS+/CD3+, INS+/CD3-,
#              INS-/CD3+) to visualize protein expression patterns and compute
#              IIRS centroids for trajectory analysis.
#
# Input: - output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds
#        - output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds
# Output: - Heatmap plots in output/RD4-islet_immune_response_signature/
#         - output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)
library(ComplexHeatmap)


# Import ------------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load batch-corrected protein abundance data
m <- readRDS("output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds")

# Load final IIRS protein set
im_sig <- readRDS("output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds")


# Process -----------------------------------------------------------------

# Require 50% data completeness per donor
m2 <- completeness_filt(m, 0.5)


# Immune Signature Heatmap ------------------------------------------------

# Subset to IIRS proteins
m_iirs <- m2[match(im_sig$UniProtAcc, fData(m2)$UniProtAcc), ]

# Split by case/control status
m_split <- split(m_iirs, "class")
m_case <- m_split[["case"]]

# Split case samples by donor
m_case$donor <- droplevels(m_case$donor)
m_case_donors <- split(m_case, "donor")

# Adjust case donors: center 25% quantile to zero
correct_case <- \(m_lp, idx) {
	f_meds <- apply(exprs(m_lp), 1, \(x) {
		quantile(x, 0.25, na.rm = T)
	})
	adjust_factor <- median(f_meds)
	exprs(m_lp) <- sweep(exprs(m_lp), 1, adjust_factor, FUN = "-")
	return(m_lp)
}

m_case_donors_corrected <- imap(m_case_donors, correct_case)


# Recombine case donors with controls
combine_m <- \(m1, m2) {
	x <- MSnbase::combine(m1, m2)
	attr(pData(x)$date, "tzone") <- attr(pData(m)$date, "tzone")
	attr(pData(x)$ship_date, "tzone") <- attr(pData(m)$ship_date, "tzone")
	return(x)
}

m_case_corrected <- reduce(m_case_donors_corrected, combine_m)
m_iirs <- suppressWarnings(
	combine_m(m_case_corrected, m_split[["control"]])
)
m_iirs$donor <- fct_relevel(m$donor, levels(m$donor))

# Calculate IIRS centroid per sample
pData(m_iirs)$im_sig_centroid <- apply(exprs(m_iirs), 2, \(x) {
	median(x, na.rm = TRUE)
})

saveRDS(
	m_iirs,
	"output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds"
)


# Center features for visualization
exprs(m_iirs) <- t(scale(t(exprs(m_iirs)), center = T, scale = F))


# Heatmap options
ht_opt$COLUMN_ANNO_PADDING = unit(3, "mm")
ht_opt$merge_legends = T
ht_opt$message = FALSE

# Create heatmap for each donor
donors <- unique(pData(m_iirs)$donor)
h_hor <- vector("list", length(donors))
names(h_hor) <- donors

# Color scales
expr_col <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
cd3_colors <- c("Positive" = "red", "Negative" = "blue")

for (i in donors) {
	d1 <- m_iirs[, pData(m_iirs)$donor == i]
	pData(d1)$CD3 <- factor(pData(d1)$CD3, levels = c("Positive", "Negative"))
	IIRS <- scale(pData(d1)$im_sig_centroid, center = T, scale = F)

	# Create annotations
	show_names <- i == "6539"
	ca1 <- HeatmapAnnotation(
		CD3 = pData(d1)$CD3,
		col = list(CD3 = cd3_colors),
		show_annotation_name = show_names,
		simple_anno_size = unit(3.75, "mm")
	)
	ca2 <- HeatmapAnnotation(
		IIRS = IIRS,
		col = list(IIRS = expr_col),
		show_annotation_name = show_names,
		simple_anno_size = unit(3.75, "mm")
	)

	# Create heatmap with conditional legend
	h_hor[[i]] <- Heatmap(
		exprs(d1),
		top_annotation = c(ca1, ca2),
		height = nrow(exprs(m_iirs)) * unit(4, "mm"),
		clustering_distance_columns = "euclidean",
		clustering_method_columns = "average",
		cluster_rows = F,
		show_column_names = F,
		column_title = i,
		show_heatmap_legend = i == "6539",
		heatmap_legend_param = if (i == "6539") list(title = "Abundance") else NULL,
		col = expr_col
	)
}

png(
	filename = "output/RD4-islet_immune_response_signature/RD4c_1-islet_immune_response_signature_heatmap.png",
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
