# ==============================================================================
# Immune Signature Differential Expression Analysis
# ==============================================================================
# Script: 5c_1-immune_signature_differential_analysis.R
# Description: Performs differential protein expression analysis associating
#              protein abundance with immune signature centroid levels using limma.
#              Filters low insulin (INS < -2) islets, applies 50% completeness
#              filtering.
#
# Input: - output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds
# Output: - output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv
# ==============================================================================

library(ComplexHeatmap)
library(MSnSet.utils)
library(tidyverse)


# Import -----------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load protein abundance with IIRS centroid
m <- readRDS("output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds")


# Process ----------------------------------------------------------------

# QC: visualize GCG and INS intensities by donor
boxplot(exprs(m)["GCG", ] ~ m$donor, xlab = "donor", ylab = "GCG")
boxplot(exprs(m)["INS", ] ~ m$donor, xlab = "donor", ylab = "INS")

# Check INS-low islet counts per donor
table(m$donor, exprs(m)["INS", ] < -2)


# Remove INS-low islets (likely non-beta cell enriched)
keep <- exprs(m)["INS", ] >= -2
m2 <- m[, keep]

# Require 50% completeness per donor
m2 <- completeness_filt(m2, 0.5)


# Z-score IIRS centroid within donors for comparable effect sizes
pData(m2) <- pData(m2) %>%
	mutate(.by = donor, im_zscore = scale(im_sig_centroid)[, 1L])


## Differential Analysis ----
res <- lapply(levels(m2$donor), function(donor_i) {
	limma_gen(m2[, m2$donor == donor_i], "~ im_zscore", "im_zscore") %>%
		tibble::rownames_to_column("protein")
}) %>%
	setNames(levels(m2$donor)) %>%
	bind_rows(.id = "donor") %>%
	mutate(donor = factor(donor, levels = unique(donor))) %>%
	select(-B, -starts_with("SE", ignore.case = FALSE)) %>%
	mutate(ZScore = qnorm(P.Value / 2, lower.tail = FALSE) * sign(logFC))

table(res$donor, res$adj.P.Val < 0.05)

# Save
res %>%
	pivot_wider(
		id_cols = protein,
		names_from = donor,
		values_from = c(logFC, ZScore, P.Value, adj.P.Val)
	) %>%
	write.csv(
		file = "output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv",
		quote = FALSE,
		row.names = FALSE
	)
