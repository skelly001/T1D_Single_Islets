# ==============================================================================
# Gene Ontology Term Heatmaps for Immune Signature Analysis
# ==============================================================================
# Script: 5f_1-GO_term_heatmaps.R
# Description: Generates detailed heatmaps for specific Gene Ontology biological
#              processes of interest showing donor-specific protein expression
#              patterns. Creates four heatmaps for: (1) extracellular matrix
#              organization, (2) hyaluronan metabolic process, (3) glycosaminoglycan
#              catabolic process, and (4) interleukin-10 production. Features are
#              selected based on significant differences between mAAb+ and
#              control donors (t-test p < 0.05).
#
# Input: - output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv
#        - data/GOBP_gene_sets.rds
# Output: - output/RD5-islet_immune_response_signature_analysis/RD5f_1-ECM-related_GO_term_heatmap.png
#         - output/RD5-islet_immune_response_signature_analysis/RD5f_1-hyaluronan_metabolic_process_GO_term_heatmap.png
#         - output/RD5-islet_immune_response_signature_analysis/RD5f_1-glycosaminoglycan_catabolic_process_GO_term_heatmap.png
#         - output/RD5-islet_immune_response_signature_analysis/RD5f_1-interleukin-10_production_GO_term_heatmap.png
# ==============================================================================

library(TMSig)
library(ComplexHeatmap)
library(tidyverse)


# Import -----------------------------------------------------------------

# Load limma differential expression results
limma_res <- read.csv(
	"output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv"
)

# Load GO Biological Process gene sets
gene_sets <- readRDS("data/GOBP_gene_sets.rds")


# Process ----------------------------------------------------------------

x <- limma_res %>%
	select(protein, starts_with("logFC")) %>%
	tibble::column_to_rownames("protein") %>%
	as.matrix()

colnames(x) <- sub("logFC_", "", colnames(x))

# Filter proteins with significant T1D vs control difference ----
group <- factor(rep(c("g1", "g2"), each = 3L))

# T-test filtering (p < 0.05)
keep <- apply(x, 1L, function(xi) {
	t.test(xi ~ group)$p.value
})

keep <- keep < 0.05
x2 <- x[keep, ]


# Clean gene set names
names(gene_sets) <- gsub(",", "", names(gene_sets))

gene_sets_filtered <- filterSets(
	x = gene_sets,
	background = rownames(x),
	min_size = 1L
)

## Heatmap Helper Functions ----
create_heatmap <- function(mat, genes, filename) {
	mat <- mat[rownames(mat) %in% genes, ]

	r <- TMSig::extendRangeNum(
		x = mat,
		nearest = 0.1
	)

	if (all(r >= 0)) {
		breaks <- at <- c(0, r[2])
		colors <- c("white", "red")
	} else if (all(r <= 0)) {
		breaks <- at <- c(r[1], 0)
		colors <- c("blue", "white")
	} else {
		breaks <- max(abs(r)) * c(-1, 0, 1)
		at <- c(r[1], 0, r[2])
		colors <- c("blue", "white", "red")
	}

	on.exit(dev.off())

	width <- convertUnit(
		unit(240, "pt"),
		"in"
	)

	width <- as.numeric(width)

	height <- convertUnit(
		unit(14, "pt") * nrow(mat),
		"in"
	)

	height <- as.numeric(height) + 1

	height <- max(height, 3.5)

	png(filename = filename, width = width, height = height, units = "in", res = 600)

	ht <- Heatmap(
		matrix = mat,
		col = circlize::colorRamp2(
			breaks = breaks,
			colors = colors
		),
		border = TRUE,
		cluster_columns = FALSE,
		column_names_side = "top",
		width = ncol(mat) * unit(14, "pt"),
		height = nrow(mat) * unit(14, "pt"),
		heatmap_legend_param = list(
			title = "logFC",
			at = at
		)
	)

	draw(
		object = ht,
		align_heatmap_legend = "heatmap_top"
	)
}

get_genes <- function(set_names, background = rownames(x)) {
	gene_sets[set_names] %>%
		unlist(use.names = FALSE) %>%
		intersect(background)
}

## ECM Organization Heatmap ----
genes1 <- get_genes(
	c(
		"extracellular matrix organization",
		"integrin-mediated signaling pathway",
		"cell-matrix adhesion"
	),
	background = rownames(x)
)


create_heatmap(
	mat = x2,
	genes = genes1,
	filename = "output/RD5-islet_immune_response_signature_analysis/RD5f_1-ECM-related_GO_term_heatmap.png"
)


## Hyaluronan Metabolism Heatmap ----
genes2 <- get_genes(
	"hyaluronan metabolic process",
	background = rownames(x)
)

create_heatmap(
	mat = x,
	genes = genes2,
	filename = "output/RD5-islet_immune_response_signature_analysis/RD5f_1-hyaluronan_metabolic_process_GO_term_heatmap.png"
)


## Glycosaminoglycan Catabolism Heatmap ----
genes3 <- get_genes(
	"glycosaminoglycan catabolic process",
	background = rownames(x)
)

create_heatmap(
	mat = x,
	genes = genes3,
	filename = "output/RD5-islet_immune_response_signature_analysis/RD5f_1-glycosaminoglycan_catabolic_process_GO_term_heatmap.png"
)


## IL-10 Production Heatmap ----
genes4 <- get_genes(
	"interleukin-10 production",
	background = rownames(x)
)

create_heatmap(
	mat = x,
	genes = genes4,
	filename = "output/RD5-islet_immune_response_signature_analysis/RD5f_1-interleukin-10_production_GO_term_heatmap.png"
)
