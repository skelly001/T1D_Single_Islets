# ==============================================================================
# WGCNA Network Analysis - Donor 6450 (Stage 1 T1D)
# ==============================================================================
# Script: 3a_1-WGCNA_donor_6450.R
# Description: Performs Weighted Gene Co-expression Network Analysis (WGCNA) on
#              proteomics data from Stage 1 T1D donor 6450. Constructs protein
#              co-expression networks, identifies modules using hierarchical
#              clustering, merges similar modules, performs over-representation
#              analysis (ORA) for GO biological processes, and correlates module
#              eigengenes with phenotypic traits (INS intensity, CD3 status).
#
# Input: - output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds
# Output: - output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_module_membership.xlsx
#         - output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_ORA_BP_results.xlsx
#         - QC and visualization plots in output/RD3-WGCNA/donor_6450/
# ==============================================================================

library(writexl)
library(WGCNA)
options(stringsAsFactors = FALSE)
disableWGCNAThreads()
library(RColorBrewer)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(tidyverse)
library(MSnSet.utils)
library(clusterProfiler)


# Import -----------------------------------------------------------------

source("data/GSEA_functions.R")

# Load protein abundance data
m <- readRDS("output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds")


# Process ----------------------------------------------------------------

m_split <- split(m, "donor")
m3 <- m_split[["6450"]]


# Remove features without gene symbols
m3 <- m3[!is.na(fData(m3)$gene_name), ]


# Require 50% data completeness
m3 <- filter_by_occurrence(m3, 0.5)

# Transpose expression matrix for WGCNA
datExpr <- t(MSnbase::exprs(m3))

## Network Construction and Module Detection ----

### Adjacency Matrix ----
# Soft power selection using biweight midcorrelation

fs::dir_create("output/RD3-WGCNA/donor_6450")

if (!file.exists("output/RD3-WGCNA/donor_6450/SoftPower_plot.png")) {
	powers <- 5:30
	sft <- pickSoftThreshold(
		data = datExpr,
		powerVector = powers,
		networkType = "signed",
		corFnc = WGCNA::bicor,
		corOptions = list(use = "pairwise.complete.obs"),
		verbose = 0,
		RsquaredCut = 0.90
	)

	png(
		str_glue("output/RD3-WGCNA/donor_6450/SoftPower_plot.png"),
		width = 10,
		height = 5,
		units = "in",
		res = 200
	)

	par(mfrow = c(1, 2))
	cex1 = 0.9
	plot(
		sft$fitIndices[, 1],
		-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
		xlab = "Soft Threshold (power)",
		ylab = "Scale Free Topology Model Fit,signed R^2",
		type = "n",
		main = paste("Scale independence")
	)
	text(
		sft$fitIndices[, 1],
		-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
		labels = powers,
		cex = cex1,
		col = "red"
	)
	# R-squared threshold line
	abline(h = 0.90, col = "red")
	plot(
		sft$fitIndices[, 1],
		sft$fitIndices[, 5],
		xlab = "Soft Threshold (power)",
		ylab = "Mean Connectivity",
		type = "n",
		main = paste("Mean connectivity")
	)
	text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
	dev.off()
}


softPower <- 15


# Compute adjacency matrix with signed network
adjacency.mat <- adjacency(
	datExpr = datExpr,
	power = softPower,
	type = "signed",
	corFnc = WGCNA::bicor,
	corOptions = list(use = "pairwise.complete.obs")
)


### Topological Overlap Matrix and Dissimilarity ----

# Compute TOM and dissimilarity
TOM <- TOMsimilarity(adjacency.mat, TOMType = "unsigned")
dissTOM <- (1 - TOM)


### Hierarchical Clustering ----

# Cluster proteins by TOM dissimilarity
geneTree <- hclust(as.dist(dissTOM), method = "average")

plot(
	geneTree,
	xlab = "",
	sub = "",
	ylab = "Dissimilarity",
	main = "Protein clustering on TOM-based dissimilarity",
	labels = FALSE,
	hang = 0.04
)


# Dynamic tree cut with minimum module size of 20
minModuleSize <- 20
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, minClusterSize = minModuleSize)


# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods, colorSeq = standardColors())


# Plot the dendrogram and colors underneath
plotDendroAndColors(
	geneTree,
	dynamicColors,
	"Dynamic Tree Cut",
	dendroLabels = FALSE,
	hang = 0.03,
	addGuide = TRUE,
	guideHang = 0.05,
	main = "Protein dendrogram and module colors"
)


### Module Merging ----

# Calculate module eigengenes
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes


# Calculate eigengene dissimilarity and cluster
MEDiss <- (1 - cor(MEs))
METree <- hclust(as.dist(MEDiss), method = "average")

# Module merge threshold: 0.2 dissimilarity (0.8 correlation)
MEDissThres <- 0.2

plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = MEDissThres, col = "red")

# Merge modules with correlation >= 0.8
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)


moduleColors <- merge$colors


# Combine feature data with module assignments
df.final <- cbind(fData(m3), moduleColors, MSnbase::exprs(m3))


temp <- sampleNames(m3)
temp <- paste0(temp, "_INS_", m3$`INS`, "_CD3_", m3$`CD3`)
names(temp) <- sampleNames(m3)

colnames(df.final)[match(sampleNames(m3), colnames(df.final))] <- temp[sampleNames(m3)]

df.final %>%
	select(-all_of(temp)) %>%
	write_xlsx(path = glue::glue("output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_module_membership.xlsx"))


## Over-Representation Analysis ----

# Map gene symbols to Entrez IDs per module
gene_clusters <- split(df.final$gene_name, df.final$moduleColors) %>%
	lapply(intersect, mappedRkeys(org.Hs.egSYMBOL2EG)) %>%
	lapply(function(mod_i) {
		sapply(as.list(org.Hs.egSYMBOL2EG[mod_i]), `[`, 1)
	}) %>%
	.[names(.) != "grey"] # Remove grey module

# Perform GO enrichment analysis for each module
ora_results <- compareCluster(
	geneClusters = gene_clusters,
	fun = "enrichGO",
	ont = "BP",
	OrgDb = "org.Hs.eg.db",
	pvalueCutoff = 1,
	qvalueCutoff = 1,
	pAdjustMethod = "BH",
	readable = TRUE,
	universe = unlist(gene_clusters)
)

# Export complete enrichment results
ora_results@compareClusterResult %>%
	filter(p.adjust < 0.05) %>%
	arrange(Cluster, p.adjust) %>%
	write_xlsx(path = str_glue("output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_ORA_BP_results.xlsx"))

# Filter significant terms
res <- ora_results@compareClusterResult %>%
	filter(p.adjust < 0.05) %>%
	group_by(Cluster) %>%
	arrange(p.adjust)

# Select top 5 terms per module for visualization
top_terms <- res %>%
	group_by(Cluster) %>%
	slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
	pull(ID) %>%
	unique()

# Format for plotting
res <- res %>%
	filter(ID %in% top_terms) %>%
	format_description(70) %>%
	arrange(p.adjust) %>%
	mutate(
		Description = factor(Description, levels = unique(Description)),
		Cluster = factor(Cluster, levels = unique(Cluster)),
		p.adjust = ifelse(p.adjust <= 1e-06, 1e-06, p.adjust) # Floor p-values for visualization
	)

# Create enrichment dotplot
p1 <- ggplot(res) +
	geom_point(aes(x = Cluster, y = Description, color = p.adjust), size = 2.7) +
	scale_x_discrete(name = "Module") +
	scale_y_discrete(name = NULL, limits = rev) +
	scale_color_viridis_c(
		name = "Adjusted\np-value",
		trans = "log10",
		breaks = trans_breaks("log10", function(x) 10^x, n = 5),
		labels = label_scientific()
	) +
	theme_bw(base_size = 8) +
	theme(
		text = element_text(color = "black", size = 8),
		axis.text.y = element_text(color = "black", size = 8),
		axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 8)
	)

# Export dotplot
ggsave(
	"output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_ORA_dotplot.png",
	p1,
	height = 8.6,
	width = 7.15,
	dpi = 600
)


### Module-Trait Correlation ----

# Add INS protein intensity as phenotype
pData(m3)$INS_intensity <- MSnbase::exprs(m3)[fData(m3)$gene_name == "INS", ]

# Create MSnSet from module eigengenes
ME_final <- merge$newMEs[, colnames(merge$newMEs) != "MEgrey"]
ME_final <- ME_final[sampleNames(m3), ]
m_me <- MSnSet(exprs = t(as.matrix(ME_final)), pData = pData(m3))
featureNames(m_me) <- str_replace(featureNames(m_me), "ME", "")

# Test module-trait associations using linear models
ptype_to_test <- c("CD3", "INS_intensity")
lm_res <- vector("list", length(ptype_to_test))
names(lm_res) <- ptype_to_test

for (i in ptype_to_test) {
	lm_res[[i]] <- limma_gen(m_me, paste0("~ ", i), i) %>%
		select(logFC, P.Value) %>%
		rownames_to_column()
}

# Combine results into data frame
lm_df <- lm_res %>%
	enframe() %>%
	unnest(value)

# Create correlation matrix from logFC values
cor_mat <- lm_df %>%
	select(-P.Value) %>%
	pivot_wider(names_from = name, values_from = logFC) %>%
	column_to_rownames() %>%
	as.matrix()

# Create adjusted p-value matrix
p_mat <- lm_df %>%
	select(-logFC) %>%
	pivot_wider(names_from = name, values_from = P.Value) %>%
	column_to_rownames() %>%
	as.matrix() %>%
	p.adjust(method = "BH") %>%
	matrix(dim(cor_mat)[[1]], dim(cor_mat)[[2]])
colnames(p_mat) <- colnames(cor_mat)
rownames(p_mat) <- rownames(cor_mat)

# Convert p-values to significance stars for heatmap annotation
comb_signif <- gtools::stars.pval(p_mat)

# Create color mapping for module annotations
anno_colors <- list(module = rownames(cor_mat))
names(anno_colors$module) <- rownames(cor_mat)

# Build custom legend for significance levels
gb <- textGrob(
	paste(".   < 0.1", "*   < 0.05", "**  < 0.01", "*** < 0.001", sep = "\n"),
	x = 0,
	y = 0,
	gp = gpar(fontsize = 6),
	hjust = 0,
	vjust = 0
)

lgd <- Legend(title = "Adjusted\np-value", grob = gb, title_gp = gpar(fontsize = 7))

# Calculate module sizes for annotation
mod_count <- table(moduleColors)
mods <- names(mod_count)

# Calculate optimal text color (black/white) for each module color background
pt_color <- sapply(rownames(cor_mat), function(mod_i) {
	c("black", "white")[
		1 + (sum(col2rgb(mod_i) * c(299, 587, 114)) / 1000 < 123)
	]
})

# Create row annotation showing module sizes
ra <- HeatmapAnnotation(
	"Module Size" = anno_simple(
		rownames(cor_mat),
		pch = as.character(mod_count[rownames(cor_mat)]),
		pt_gp = gpar(col = pt_color),
		col = anno_colors$module,
		width = unit(0.2, "in"),
		pt_size = unit(8, "points")
	),
	which = "row",
	annotation_name_gp = gpar(fontsize = 10),
	show_legend = FALSE
)

# Build module-trait correlation heatmap
ht <- Heatmap(
	cor_mat,
	cluster_columns = FALSE,
	border = TRUE,
	row_title_gp = gpar(fontsize = 10),
	column_title_gp = gpar(fontsize = 10),
	row_dend_side = "right",
	row_names_side = "left",
	column_title_side = "top",
	column_dend_height = unit(0.2, "in"),
	row_dend_width = unit(0.2, "in"),
	col = circlize::colorRamp2(breaks = c(-0.2, 0, 0.2), colors = c("blue", "white", "red")),
	column_names_gp = gpar(fontsize = 10),
	row_names_gp = gpar(fontsize = 10),
	heatmap_legend_param = list(
		title = "logFC",
		title_gp = gpar(fontsize = 8),
		title_position = "lefttop-rot",
		legend_height = unit(0.85, "in"),
		grid_width = unit(0.1, "in"),
		labels_gp = gpar(fontsize = 8)
	),
	left_annotation = ra,
	cell_fun = function(j, i, x, y, width, height, fill) {
		grid.text(comb_signif[i, j], x, y, gp = gpar(fontsize = 8))
	}
)

# Export heatmap
png(
	str_glue("output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_module_correlation_heatmap.png"),
	width = 2.5,
	height = 3.6,
	units = "in",
	res = 600
)
draw(ht, heatmap_legend_list = lgd, merge_legends = TRUE)
dev.off()
