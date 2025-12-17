# ==============================================================================
# Protein Co-Expression Cluster Heatmap Visualization
# ==============================================================================
# Script: 8c_1-cluster_heatmaps.R
# Description: Creates heatmaps for protein co-expression clusters with Gene Ontology
#              annotations. Calculates module eigengenes, correlates clusters with
#              IIRS and BCP centroids, and generates per-cluster heatmaps with
#              hierarchical clustering. Large clusters are split into sub-clusters,
#              and results are combined into a single PDF.
#
# Input: - output/RD8-clustering/RD8b_1-clustering_results.RData
#        - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
#        - output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds
#        - output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx
# Output: - output/RD8-clustering/RD8c_1-cluster_heatmaps.pdf
#         - output/RD8-clustering/RD8c_1-LGALS3_LGALS3BP_heatmap.pdf
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)
library(ComplexHeatmap)


# Import ------------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load clustering results
(load("output/RD8-clustering/RD8b_1-clustering_results.RData"))

# Load protein abundance data with centroids
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")

# Load IIRS data
m_imc <- readRDS(
	"output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds"
)

# Load BCP protein features
bcp_proteins <- openxlsx::read.xlsx(
	"output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx"
)


# Process -----------------------------------------------------------------

# Subset to genes used in clustering
m2 <- m[fData(m)$gene_name %in% WGCNA_Summary[[4]]$Gene, ]

# Add cluster assignments to feature data
clust_assignment <- WGCNA_Summary[[4]][, c("Gene", "Cluster")]
fData(m2) <- fData(m2) %>%
	rownames_to_column() %>%
	left_join(clust_assignment, by = c("gene_name" = "Gene")) %>%
	column_to_rownames()


## Cluster-Trait Correlation ----

# Calculate module eigengenes (case donors only)
exprDat <- t(exprs(m2)[, pData(m2)$class == "case"])
MEList <- WGCNA::moduleEigengenes(exprDat, colors = fData(m2)$Cluster)
MEs <- MEList$eigengenes


# Prepare IIRS and BCP for correlation
case_donor_ids <- sampleNames(m2)[pData(m2)$class == "case"]

IIRS <- pData(m_imc) %>%
	rownames_to_column() %>%
	mutate(im_sig_centroid = scale(im_sig_centroid, center = T, scale = F)) %>%
	pull(im_sig_centroid, rowname)
IIRS_case <- IIRS[case_donor_ids]

m_bcp <- m2[featureNames(m2) %in% bcp_proteins$gene_name, ]
BCP <- apply(exprs(m_bcp), 2, \(x) {
	median(x, na.rm = T)
})
BCP <- scale(BCP, center = T, scale = F)[, 1]
BCP_case <- BCP[case_donor_ids]


# Correlate clusters with IIRS
IIRS_c <- apply(MEs, 2, \(x) {
	cor.test(IIRS_case, x, use = "pairwise.complete.obs")$estimate
})
IIRS_p <- apply(MEs, 2, \(x) {
	cor.test(IIRS_case, x, use = "pairwise.complete.obs")$p.value
})
IIRS_cor <- tibble(
	cluster = as.integer(str_extract(names(IIRS_c), "\\d+")),
	cor = IIRS_c,
	pval = IIRS_p
)

# Correlate clusters with BCP
BCP_c <- apply(MEs, 2, \(x) {
	cor.test(BCP_case, x, use = "pairwise.complete.obs")$estimate
})
BCP_p <- apply(MEs, 2, \(x) {
	cor.test(BCP_case, x, use = "pairwise.complete.obs")$p.value
})
BCP_cor <- tibble(
	cluster = as.integer(str_extract(names(BCP_c), "\\d+")),
	cor = BCP_c,
	pval = BCP_p
)

# Center features
exprs(m2) <- t(scale(t(exprs(m2)), center = T, scale = F))


## Generate Cluster Heatmaps ----------------------------------------------

ht_opt$COLUMN_ANNO_PADDING = unit(3, "mm")
ht_opt$merge_legends = T
ht_opt$message = FALSE


# Create temporary output directory
dir.create("output/RD8-clustering/temp")

# Define donors for iteration
donors <- unique(pData(m2)$donor)


# Identify large clusters requiring subdivision
clust_too_long <- fData(m2) %>%
	group_by(Cluster) %>%
	summarise(n = n()) %>%
	arrange(desc(n)) %>%
	filter(n > 1000)

m_to_split <- m2[fData(m2)$Cluster %in% clust_too_long$Cluster, ]
clust48_dist <- dist(exprs(m_to_split)[, pData(m_to_split)$donor == "6450"])
dend <- hclust(clust48_dist)
split_clust <- dendextend::cutree(dend, h = 15)
split_clust <- ifelse(split_clust == 1, "A", "B") %>%
	enframe("gene_name", "sub_clust")

# Annotate sub-cluster assignments
fData(m2) <- fData(m2) %>%
	rownames_to_column() %>%
	left_join(split_clust, by = "gene_name") %>%
	mutate(Cluster = ifelse(!is.na(sub_clust), paste0(Cluster, sub_clust), Cluster)) %>%
	select(-sub_clust) %>%
	column_to_rownames()


# Generate heatmaps for all clusters
clusters <- sort(unique(fData(m2)$Cluster))
expr_col <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

for (cluster in clusters) {
	m_lp <- m2[fData(m2)$Cluster == cluster, ]

	# horizontal heatmap groups
	h_hor <- vector("list", 6)
	names(h_hor) <- donors

	for (i in donors) {
		donor_ids <- sampleNames(m2)[m2$donor == i]
		d1 <- m_lp[, donor_ids]
		featureNames(d1) <- str_extract(featureNames(d1), "^[^\\.]+")

		# annotations
		d_IIRS <- IIRS[donor_ids]
		d_BCP <- BCP[donor_ids]

		# Create annotations
		show_names <- i == "6539"
		ca1 <- HeatmapAnnotation(
			IIRS = d_IIRS,
			col = list(IIRS = expr_col),
			show_annotation_name = show_names,
			simple_anno_size = unit(3.75, "mm")
		)

		ca2 <- HeatmapAnnotation(
			BCP = d_BCP,
			col = list(BCP = expr_col),
			show_annotation_name = show_names,
			simple_anno_size = unit(3.75, "mm")
		)

		# Determine parameters based on donor ID
		show_legend <- FALSE
		legend_param <- NULL

		if (i == "6539") {
			show_legend <- TRUE
			legend_param <- list(title = "Abundance")
		}

		h_hor[[i]] <- Heatmap(
			exprs(d1),
			top_annotation = c(ca1, ca2),
			clustering_method_columns = "average",
			clustering_method_rows = "average",
			show_column_names = F,
			height = nrow(exprs(d1)) * unit(4, "mm"),
			column_title = i,
			show_heatmap_legend = show_legend,
			heatmap_legend_param = legend_param,
			col = expr_col
		)
	}

	# hmap title
	cluster_num <- str_extract(cluster, "\\d{1,2}")
	title_components <- WGCNA_Summary[[3]][
		WGCNA_Summary[[3]]$GO_Cluster == as.character(cluster_num),
		c("p_value", "term_name")
	]
	p_val <- ifelse(
		nrow(title_components) == 1,
		as.character(format(title_components$p_value, digit = 2, scientific = T)),
		"NS"
	)
	representative_GO_term <- ifelse(nrow(title_components) == 1, title_components$term_name, "NS")
	IIRS_c <- IIRS_cor[IIRS_cor$cluster == cluster_num, "cor", drop = T] %>%
		format(digit = 2, scientific = F)
	IIRS_p <- IIRS_cor[IIRS_cor$cluster == cluster_num, "pval", drop = T] %>%
		format(digit = 2, scientific = T)
	BCP_c <- BCP_cor[BCP_cor$cluster == cluster_num, "cor", drop = T] %>%
		format(digit = 2, scientific = F)
	BCP_p <- BCP_cor[BCP_cor$cluster == cluster_num, "pval", drop = T] %>%
		format(digit = 2, scientific = T)
	top_go_term <- str_wrap(
		str_glue("Top GO term: {representative_GO_term} (p = {p_val})"),
		width = 90
	)

	pdf(
		file = str_glue("output/RD8-clustering/temp/Cluster_{cluster}_pval_{p_val}.pdf"),
		onefile = T,
		height = 3 + 0.16 * nrow(exprs(m_lp)),
		width = 11
	)
	draw(
		h_hor[[1]] + h_hor[[2]] + h_hor[[3]] + h_hor[[4]] + h_hor[[5]] + h_hor[[6]],
		gap = unit(5, "mm"),
		column_title = str_glue(
			"Cluster: {cluster}
                                   {top_go_term}
                                    IIRS Cor: {IIRS_c} (p = {IIRS_p})
                                    BCP Cor: {BCP_c} (p = {BCP_p})"
		)
	)
	dev.off()
}


## Combine PDFs ----
pdfs <- list.files("./output/RD8-clustering/temp", full.names = T)

file_path_ordered <- tibble(file_path = pdfs) %>%
	mutate(num = str_extract(file_path, "(?<=Cluster_)\\d{1,2}"), num = as.numeric(num)) %>%
	arrange(num, file_path)

pdftools::pdf_combine(
	file_path_ordered$file_path,
	output = "output/RD8-clustering/RD8c_1-cluster_heatmaps.pdf"
)

# Clean up temporary files
unlink("output/RD8-clustering/temp", recursive = TRUE)


## LGALS3/LGALS3BP Heatmap ------------------------------------------------

m_lp <- m2[grepl("LGALS3", fData(m2)$gene_name), ]


# horizontal heatmap groups
h_hor <- vector("list", 6)
names(h_hor) <- donors

expr_col <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

for (i in donors) {
	donor_ids <- sampleNames(m2)[m2$donor == i]
	d1 <- m_lp[, donor_ids]
	featureNames(d1) <- str_extract(featureNames(d1), "^[^\\.]+")

	# annotations
	d_IIRS <- IIRS[donor_ids]
	d_BCP <- BCP[donor_ids]

	# Create annotations
	show_names <- i == "6539"
	ca1 <- HeatmapAnnotation(
		IIRS = d_IIRS,
		col = list(IIRS = expr_col),
		show_annotation_name = show_names,
		simple_anno_size = unit(3.75, "mm")
	)

	ca2 <- HeatmapAnnotation(
		BCP = d_BCP,
		col = list(BCP = expr_col),
		show_annotation_name = show_names,
		simple_anno_size = unit(3.75, "mm")
	)

	# Determine parameters based on donor ID
	show_legend <- FALSE
	legend_param <- NULL

	if (i == "6539") {
		show_legend <- TRUE
		legend_param <- list(title = "Abundance")
	}

	h_hor[[i]] <- Heatmap(
		exprs(d1),
		top_annotation = c(ca1, ca2),
		clustering_method_columns = "average",
		clustering_method_rows = "average",
		height = nrow(exprs(m_lp)) * unit(4, "mm"),
		show_column_names = F,
		column_title = i,
		show_heatmap_legend = show_legend,
		heatmap_legend_param = legend_param,
		col = expr_col
	)
}


# hmap title
IIRS_c <- IIRS_cor[IIRS_cor$cluster == cluster_num, "cor", drop = T] %>%
	format(digit = 2, scientific = F)
IIRS_p <- IIRS_cor[IIRS_cor$cluster == cluster_num, "pval", drop = T] %>%
	format(digit = 2, scientific = T)
BCP_c <- BCP_cor[BCP_cor$cluster == cluster_num, "cor", drop = T] %>%
	format(digit = 2, scientific = F)
BCP_p <- BCP_cor[BCP_cor$cluster == cluster_num, "pval", drop = T] %>%
	format(digit = 2, scientific = T)


pdf(
	file = str_glue("output/RD8-clustering/RD8c_1-LGALS3_LGALS3BP_heatmap.pdf"),
	onefile = T,
	height = 3 + 0.16 * nrow(exprs(m_lp)),
	width = 11
)
draw(
	h_hor[[1]] + h_hor[[2]] + h_hor[[3]] + h_hor[[4]] + h_hor[[5]] + h_hor[[6]],
	gap = unit(5, "mm"),
	column_title = str_glue(
		"LGALS3/LGALS3BP
                                    IIRS Cor: {IIRS_c} (p = {IIRS_p})
                                    BCP Cor: {BCP_c} (p = {BCP_p})"
	)
)
dev.off()
