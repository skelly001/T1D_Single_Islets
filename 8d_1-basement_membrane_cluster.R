# ==============================================================================
# Basement Membrane Cluster Correlation Analysis
# ==============================================================================
# Script: 8d_1-basement_membrane_cluster.R
# Description: Correlates the basement membrane co-expression cluster (Cluster 67)
#              with beta cell profile (BCP) and islet immune response signature
#              (IIRS) centroids. Loads WGCNA clustering results and protein
#              abundance data, computes module eigengenes, removes islet selection
#              outliers, and plots the Cluster 67 eigengene against
#              BCP and IIRS per donor with correlation statistics.
#
# Input: - output/RD8-clustering/RD8b_1-clustering_results.RData
#        - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
# Output: - output/RD8-clustering/RD8d_1-basement_membrane_vs_BCP.png
#         - output/RD8-clustering/RD8d_1-basement_membrane_vs_IIRS.png
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)


# Import -----------------------------------------------------------------

# Clustered proteins
(load("output/RD8-clustering/RD8b_1-clustering_results.RData"))

# Protein abundance data with centroids
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")


# Process -----------------------------------------------------------------

# Output location
base_path <- "output/RD8-clustering"
fs::dir_create(base_path)

# Subset to genes used in clustering
m2 <- m[fData(m)$gene_name %in% WGCNA_Summary[[4]]$Gene, ]

# Add cluster assignments to feature data
clust_assignment <- WGCNA_Summary[[4]][, c("Gene", "Cluster")]
fData(m2) <- fData(m2) %>%
	rownames_to_column() %>%
	left_join(clust_assignment, by = c("gene_name" = "Gene")) %>%
	column_to_rownames()


# Evaluate islet selection outliers
ep_data <- as.data.frame.MSnSet(m2) %>%
	group_by(donor) %>%
	mutate(BCP_GCG_low = (bcp_centroid < -1) & (GCG < -2))


## Basement membrane proteins (ECM) vs BCP/IIRS Correlation ----

# Remove COL4A1 from clust 67 (outlier)
m2 <- m2[fData(m2)$gene_name != "COL4A1", ]
exprDat <- t(exprs(m2))

# Calculate module eigengenes
MEList <- WGCNA::moduleEigengenes(exprDat, colors = fData(m2)$Cluster)
MEs <- MEList$eigengenes

# Remove outliers
plot_data <- ep_data %>%
	left_join(rownames_to_column(MEs), by = c("sample_name" = "rowname")) %>%
	ungroup() %>%
	filter(!BCP_GCG_low)

# Plot specification
y_max <- max(plot_data$ME67, na.rm = T) + 0.04
y_min <- min(plot_data$ME67, na.rm = T)

plot_bm <- \(centroid) {
	label <- if (centroid == "bcp_centroid") "BCP" else "IIRS"
	ggplot(plot_data, aes(.data[[centroid]], ME67)) +
		facet_wrap("donor") +
		geom_point() +
		geom_smooth(method = "lm") +
		theme_bw(12) +
		ggpubr::stat_cor() +
		scale_y_continuous(limits = c(y_min, y_max)) +
		labs(x = label, y = "Basement Membrane (Cluster 67)")

	ggsave(
		paste0(base_path, "/RD8d_1-basement_membrane_vs_", label, ".png"),
		width = 8,
		height = 6,
		dpi = 600
	)
}

# Save plots
walk(c("bcp_centroid", "im_sig_centroid"), plot_bm)
