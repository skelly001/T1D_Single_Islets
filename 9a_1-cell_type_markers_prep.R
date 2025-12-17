# ==============================================================================
# Cell Type Marker Preparation from scRNA-seq Reference Data
# ==============================================================================
# Script: 9a_1-cell_type_markers_prep.R
# Description: Prepares scRNA-seq reference data from Azimuth human pancreas
#              dataset. Matches mass spectrometry-identified proteins to scRNA-seq
#              genes using UniProt database. Calculates normalized expression count
#              ratios by cell type for identifying cell type markers.
#
# Input: - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
#        - data/azimuth_human_pancreas_scRNAseq_fullref.rds
#        - data/uniprotkb_Human_2023_10_25.tsv.gz
# Output: - output/RD9-cell_type_marker_QC/RD9a_1-scRNA-Seq_count_ratios.rds
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)
library(Seurat)


# Import Data -------------------------------------------------------------

# Load protein abundance data
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")

# Load scRNA-seq reference data (Azimuth human pancreas)
fullref <- readRDS("azimuth-references/human_pancreas_snakemake/seurat_objects/fullref.Rds")

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

# Process -----------------------------------------------------------------

dir.create("output/RD9-cell_type_marker_QC")


# Extract scRNA-seq cell type metadata
meta <- fullref@meta.data %>%
	mutate(cell = rownames(.)) %>%
	select(c(cell, celltype)) %>%
	`rownames<-`(NULL)

# Observed MS proteins for scRNA-seq gene filtering
ms_proteins <- fData(m)$UniProtAcc

## Primary Gene Name Matching ----
uniprot_sp_pirimary <- uniprot %>%
	filter(Reviewed == "reviewed", !is.na(Gene.Names..primary.)) %>%
	select(Entry, Gene.Names..primary.) %>%
	filter(Entry %in% ms_proteins) %>%
	separate_longer_delim(cols = Gene.Names..primary., delim = "; ")

# First pass: match using primary gene names
sct_counts1 <- as.matrix(fullref@assays[["SCT"]]@counts) %>%
	as.data.frame() %>%
	mutate(gene = rownames(.)) %>%
	`rownames<-`(NULL) %>%
	right_join(uniprot_sp_pirimary, by = c("gene" = "Gene.Names..primary.")) %>%
	relocate(c(Entry, gene), .before = colnames(.)[1]) %>%
	mutate(missing = if_else(is.na(rowSums(.[-c(1:2)])), TRUE, FALSE), .after = gene) %>%
	arrange(missing) %>%
	distinct(Entry, .keep_all = T)

# Successful primary matches
sct_counts1_success <- sct_counts1 %>%
	filter(missing == FALSE) %>%
	select(-missing)

# Remaining unmatched proteins
sct_counts1_remain <- sct_counts1 %>%
	filter(missing == TRUE) %>%
	select(Entry)

## Alternate Gene Name Matching ----
uniprot_sp_alt <- uniprot %>%
	filter(Reviewed == "reviewed") %>%
	select(-c(Reviewed, Entry.Name, Protein.names, Organism, Length)) %>%
	filter(Entry %in% sct_counts1_remain$Entry) %>%
	unite(
		c(Gene.Names, Gene.Names..ordered.locus., Gene.Names..ORF., Gene.Names..synonym.),
		col = alt_gene_name,
		sep = "; ",
		na.rm = T
	) %>%
	separate_longer_delim(cols = alt_gene_name, delim = "; ") %>%
	separate_longer_delim(cols = alt_gene_name, delim = " ") %>%
	select(-Gene.Names..primary.) %>%
	distinct(Entry, alt_gene_name) %>%
	mutate(alt_gene_name = sub(";", "", alt_gene_name))


# Second pass: match remaining proteins using alternate gene names
sct_counts2 <- as.matrix(fullref@assays[["SCT"]]@counts) %>%
	as.data.frame() %>%
	mutate(gene = rownames(.)) %>%
	`rownames<-`(NULL) %>%
	right_join(uniprot_sp_alt, by = c("gene" = "alt_gene_name")) %>%
	relocate(c(Entry, gene), .before = colnames(.)[1]) %>%
	mutate(missing = if_else(is.na(rowSums(.[-c(1:2)])), TRUE, FALSE), .after = gene) %>%
	arrange(missing) %>%
	distinct(Entry, .keep_all = T)

# Successful alternate matches
sct_counts2_success <- sct_counts2 %>%
	filter(missing == FALSE) %>%
	select(-missing)

# Remaining unmatched (no scRNA-seq data)
sct_counts2_remain <- sct_counts2 %>%
	filter(missing == TRUE) %>%
	select(Entry, gene)


## Combine Results ----
sct_counts <- bind_rows(sct_counts1_success, sct_counts2_success, sct_counts2_remain) %>%
	select(-gene)


## Calculate Expression Ratios ----
sct_counts <- sct_counts %>%
	t() %>%
	as.data.frame() %>%
	`colnames<-`(.[1, ]) %>%
	mutate(cell = rownames(.), .before = colnames(.)[1]) %>%
	slice(2:nrow(.)) %>%
	`rownames<-`(NULL)

sct_counts_meta <- sct_counts %>%
	left_join(meta, by = "cell") %>%
	relocate(celltype, .before = cell) %>%
	select(-cell) %>%
	mutate(across(!matches("celltype"), as.numeric))

gene_distributions <- sct_counts_meta %>%
	group_by(celltype) %>%
	summarise(across(everything(), sum)) %>%
	ungroup() %>%
	filter(!is.na(celltype))

RNA_Seq_NA_genes <- gene_distributions[, colSums(is.na(gene_distributions)) != 0] %>%
	t() %>%
	as.data.frame() %>%
	rownames_to_column() %>%
	dplyr::rename(missing_genes = rowname)

gene_distributions <- gene_distributions[, colSums(is.na(gene_distributions)) == 0]

gene_ratios <- gene_distributions %>%
	mutate(across(where(~ !is.character(.x) && sum(.x) != 0), \(x) {
		x / sum(x)
	}))


# Validate islet markers (GCG, INS, SCG2)
prot_islet <- gene_ratios %>%
	select(c("celltype", "P01275", "P01308", "P13521"))

# Validate acinar markers (REG1A, REG1B, CPA1)
prot_acinar <- gene_ratios %>%
	select(c("celltype", "P05451", "P48304", "P15085"))


# Save
saveRDS(gene_ratios, file = "output/RD9-cell_type_marker_QC/RD9a_1-scRNA-Seq_count_ratios.rds")
