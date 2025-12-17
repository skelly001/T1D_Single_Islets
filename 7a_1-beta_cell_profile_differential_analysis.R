# ==============================================================================
# Beta Cell Profile Differential Expression Analysis
# ==============================================================================
# Script: 7a_1-beta_cell_profile_differential_analysis.R
# Description: Performs differential protein expression analysis for Beta Cell
#              Profile using limma. Applies 50% completeness filtering, standardizes
#              BCP centroid values (Z-score) within donors for CAMERA-PR, and performs
#              donor-specific linear modeling to identify proteins associated with
#              beta cell function.
#
# Input: - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
# Output: - output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv
# ==============================================================================

library(ComplexHeatmap)
library(MSnSet.utils)
library(tidyverse)


# Import -----------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load protein abundance data with centroids
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")


# Process ----------------------------------------------------------------

dir.create("output/RD7-beta_cell_profile_analysis")


# Require 50% completeness per donor
m2 <- completeness_filt(m, 0.5)

# Z-score BCP centroid within donors for comparable effect sizes
pData(m2) <- pData(m2) %>%
	mutate(.by = donor, im_zscore = scale(bcp_centroid)[, 1])


## Differential analysis ----
res <- lapply(levels(m2$donor), function(donor_i) {
	limma_gen(m2[, m2$donor == donor_i], "~ im_zscore", "im_zscore") %>%
		rownames_to_column("protein")
}) %>%
	setNames(levels(m2$donor)) %>%
	bind_rows(.id = "donor") %>%
	mutate(donor = factor(donor, levels = unique(donor))) %>%
	select(-B, -starts_with("SE", ignore.case = FALSE)) %>%
	mutate(ZScore = qnorm(P.Value / 2, lower.tail = FALSE) * sign(logFC))

table(res$donor, res$adj.P.Val < 0.05)

# Save
out <- res %>%
	pivot_wider(
		id_cols = protein,
		names_from = donor,
		values_from = c(logFC, ZScore, P.Value, adj.P.Val)
	)

write.csv(
	out,
	file = "output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv",
	quote = FALSE,
	row.names = FALSE
)
