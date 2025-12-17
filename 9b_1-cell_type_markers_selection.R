# ==============================================================================
# Cell Type Marker Selection Based on Expression Specificity
# ==============================================================================
# Script: 9b_1-cell_type_markers_selection.R
# Description: Selects high-specificity cell type markers from scRNA-seq reference
#              data using expression ratio thresholds. Filters proteins with >0.75
#              expression ratio in a single cell type, indicating cell type-specific
#              markers. Creates summary statistics showing number of markers per cell
#              type for islet selection quality control validation.
#
# Input: - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
#        - output/RD9-cell_type_marker_QC/RD9a_1-scRNA-Seq_count_ratios.rds
#        - data/uniprotkb_Human_2023_10_25.tsv.gz
# Output: - output/RD9-cell_type_marker_QC/RD9b_1-scRNA-Seq_cell_type_markers.rds
# ==============================================================================

library(tidyverse)


# Import Data -------------------------------------------------------------

# Load protein abundance data
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")

# Load UniProt database for gene name mapping
up_path <- "data/uniprotkb_Human_2023_10_25.tsv.gz"
uniprot <- data.table::fread(
	file = up_path,
	data.table = F,
	stringsAsFactors = F,
	showProgress = T,
	na.strings = c("", "NA"),
	check.names = T
)

# Load scRNA-seq expression ratio data
gene_ratios <- readRDS("output/RD9-cell_type_marker_QC/RD9a_1-scRNA-Seq_count_ratios.rds")


# Process -----------------------------------------------------------------

# Create gene name lookup key
key <- uniprot %>%
	select(Entry, Gene.Names..primary.) %>%
	dplyr::rename(gene_name = Gene.Names..primary.) %>%
	drop_na(gene_name)


# Transpose and annotate with gene names
gene_ratios2 <- gene_ratios %>%
	column_to_rownames("celltype") %>%
	t() %>%
	as.data.frame() %>%
	rownames_to_column("UniProtAcc") %>%
	left_join(key, by = c("UniProtAcc" = "Entry")) %>%
	mutate(featureName = paste(gene_name, UniProtAcc, sep = "_"), .keep = "unused") %>%
	column_to_rownames("featureName")


## Select High-Specificity Markers ----
threshold <- 0.75

cell_markers <- gene_ratios2 %>%
	rownames_to_column() %>%
	pivot_longer(cols = !matches("rowname"), names_to = "Cell Type", values_to = "Ratio") %>%
	filter(Ratio > threshold) %>%
	mutate(num_markers = n(), .by = `Cell Type`)

cell_markers %>%
	distinct(`Cell Type`, num_markers)

saveRDS(cell_markers, "output/RD9-cell_type_marker_QC/RD9b_1-scRNA-Seq_cell_type_markers.rds")
