# ==============================================================================
# Beta Cell Profile Pathway Enrichment Analysis Using CAMERA-PR
# ==============================================================================
# Script: 7c_1-beta_cell_profile_CAMERA-PR.R
# Description: Performs gene set enrichment analysis using CAMERA-PR (CCorrelation
#              Adjusted Mean RAnk - Pre-Ranked) to identify biological
#              pathways and processes associated with beta cell profile variation.
#              Applies gene set size filtering, Jaccard similarity clustering
#              (cutoff: 2/3), and filters for significance (FDR < 0.05) in at least
#              one mAAb+ donor.
#
# Input: - output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv
#        - data/GOBP_gene_sets.rds
# Output: - output/RD7-beta_cell_profile_analysis/RD7c_1-gene_set_table.csv
#         - output/RD7-beta_cell_profile_analysis/RD7c_1-beta_cell_profile_CAMERA-PR_significant_results.csv
#         - output/RD7-beta_cell_profile_analysis/RD7c_1-beta_cell_profile_CAMERA-PR_bubble_plot.png
# ==============================================================================

library(TMSig)
library(tidyverse)

# Import -----------------------------------------------------------------

# Load limma differential expression results
limma_res <- read.csv(
	"output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv"
)

# Load GO Biological Process gene sets
gene_sets <- readRDS("data/GOBP_gene_sets.rds")


# Process ----------------------------------------------------------------

## Prepare Data ----
x <- limma_res %>%
	pivot_longer(
		cols = -protein,
		names_to = c(".value", "donor"),
		names_pattern = "(.*)_(.*)"
	) %>%
	mutate(
		donor = paste0("donor_", donor),
		donor = factor(donor, levels = unique(donor)),
		z = sign(logFC) * qnorm(P.Value / 2, lower.tail = FALSE)
	) %>%
	pivot_wider(
		id_cols = protein,
		names_from = "donor",
		values_from = "z"
	) %>%
	column_to_rownames("protein") %>%
	as.matrix()


## Prepare gene sets ----

# Clean gene set names
names(gene_sets) <- gsub(",", "", names(gene_sets))

gene_sets_final <- filterSets(
	x = gene_sets,
	background = rownames(x),
	min_size = 5L,
	max_size = nrow(x) - 1L
)

## Export Gene Set Membership ----
gene_set_df <- vapply(
	gene_sets_final,
	function(li) {
		paste(sort(li), collapse = ";")
	},
	character(1L)
) %>%
	enframe(name = "gene_set", value = "genes")

write.csv(
	x = gene_set_df,
	file = "output/RD7-beta_cell_profile_analysis/RD7c_1-gene_set_table.csv",
	quote = FALSE,
	row.names = FALSE
)

## Run CAMERA-PR ----
res <- cameraPR.matrix(statistic = x, index = gene_sets_final, min.size = 5L) %>%
	mutate(
		donor = sub("donor_", "", Contrast),
		donor = factor(donor, levels = unique(donor)),
		Contrast = NULL
	) %>%
	select(donor, everything())

# Cluster significant gene sets by Jaccard similarity
sets_to_cluster <- res %>%
	filter(.by = GeneSet, any(FDR[donor %in% levels(donor)[1:3]] < 0.05)) %>%
	pull(GeneSet) %>%
	unique() %>%
	gene_sets_final[.]

cluster_df <- clusterSets(
	x = sets_to_cluster,
	type = "jaccard",
	cutoff = 2 / 3
) %>%
	select(GeneSet = set, cluster)

max(cluster_df$cluster)

res <- left_join(res, cluster_df, by = "GeneSet") %>%
	relocate(cluster, .before = GeneSet)

# Save
res_out <- res %>%
	# Significant in at least one donor
	filter(.by = GeneSet, any(FDR < 0.05)) %>%
	pivot_wider(
		id_cols = c(cluster, GeneSet, NGenes),
		names_from = donor,
		values_from = c(ZScore, PValue, FDR)
	)

write.csv(
	res_out,
	file = "output/RD7-beta_cell_profile_analysis/RD7c_1-beta_cell_profile_CAMERA-PR_significant_results.csv",
	quote = FALSE,
	row.names = FALSE
)


## Bubble Heatmap ----

# Load manually curated gene set selection
file <- "data/RD7c_1-beta_cell_profile_CAMERA-PR_significant_results_selected.xlsx"
selected_sets <- readxl::read_xlsx(file) %>%
	filter(selected) %>%
	pull(GeneSet)


png(
	filename = "output/RD7-beta_cell_profile_analysis/RD7c_1-beta_cell_profile_CAMERA-PR_bubble_plot.png",
	width = 8,
	height = 6.5,
	units = "in",
	res = 600
)

res %>%
	filter(GeneSet %in% selected_sets) %>%
	enrichmap(
		n_top = Inf,
		set_column = "GeneSet",
		contrast_column = "donor",
		statistic_column = "ZScore",
		padj_column = "FDR",
		padj_legend_title = "BH Adjusted\nP-Value",
		colors = c("blue", "red"),
		heatmap_args = list(
			column_names_side = "top",
			heatmap_legend_param = list(
				title = "Z-Score"
			)
		)
	)

dev.off()
