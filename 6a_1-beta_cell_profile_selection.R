# ==============================================================================
# Beta Cell Profile Selection
# ==============================================================================
# Script: 6a_1-beta_cell_profile_selection.R
# Description: Identifies and selects Beta Cell Profile (BCP) proteins
#              based on correlation analysis with key beta cell markers INS (insulin)
#              and ENTPD3. Applies 50% completeness filtering across all donors,
#              ranks proteins by mean correlation, and selects top 40 candidates
#              plus additional markers (GAD2, PCSK1) to define the beta cell profile.
#
# Input: - output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds
#        - output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds
# Output: - output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx
#         - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)


# Import ------------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load protein abundance data with IIRS
m <- readRDS("output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds")

# Load IIRS data
m_imc <- readRDS(
	"output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds"
)


# Process -----------------------------------------------------------------

dir.create("output/RD6-beta_cell_profile")

# Require 50% completeness per donor
m2 <- completeness_filt(m, 0.5)


## Beta Cell Marker Correlation ----
cor_fun <- \(x) {
	prot <- exprs(m2)[fData(m2)$gene_name == x, ]
	prot_c <- apply(t(exprs(m2)), 2, \(x) {
		cor.test(x, prot, use = "pairwise.complete.obs")$estimate
	})
	prot_p <- apply(t(exprs(m2)), 2, \(x) {
		cor.test(x, prot, use = "pairwise.complete.obs")$p.value
	})
	prot_cor <- tibble(
		gene_name = names(prot_c),
		!!str_glue("{x}_cor") := prot_c,
		!!str_glue("{x}_pval") := prot_p
	)
}

ENTPD3 <- cor_fun("ENTPD3")
INS <- cor_fun("INS")

marker_cor <- left_join(INS, ENTPD3, "gene_name")

# Rank all proteins by mean correlation
all_markers <- marker_cor %>%
	mutate(mean_cor = mean(c(INS_cor, ENTPD3_cor)), .by = "gene_name") %>%
	arrange(desc(mean_cor))

# Select top 40 candidates
marker_candidates <- marker_cor %>%
	mutate(mean_cor = mean(c(INS_cor, ENTPD3_cor)), .by = "gene_name") %>%
	slice_max(order_by = mean_cor, n = 40)


# Finalize BCP feature list (add GAD2, PCSK1; ensure INS first)
fin_feat_gn <- c("INS", setdiff(marker_candidates$gene_name, c("INS")), "GAD2", "PCSK1")

# Reorder by final selection order
fin_feat <- all_markers[match(fin_feat_gn, all_markers$gene_name), ]

# Save
openxlsx::write.xlsx(
	fin_feat,
	"output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx"
)


## Calculate BCP Centroid ----
m_sig <- m[fData(m)$gene_name %in% fin_feat$gene_name, ]
pData(m)$bcp_centroid <- apply(exprs(m_sig), 2, \(x) {
	median(x, na.rm = TRUE)
})

saveRDS(m, "output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")
