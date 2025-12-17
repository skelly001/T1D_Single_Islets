# ==============================================================================
# Preprocessing for Stage 1 T1D Donors (1, 2, 3) - iBAQ Quantification
# ==============================================================================
# Script: 2c_1-process_msnsets_donors_1_2_3_iBAQ.R
# Description: Preprocesses iBAQ (intensity-based absolute quantification) data
#              for Stage 1 T1D donors 6450, 6521, and 6267. Calculates protein-level
#              iBAQ values, performs log2 transformation, normalization, outlier
#              removal, and tissue block batch effect correction.
#
# Input: - output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData
#        - data/uniprot_database_observed_proteins.tsv
#        - data/Results_7_combined/2024-10-09-decoys-contam-uniprotkb_Human_AND_reviewed_true_AND_m_2024_10_09.fasta.fas
# Output: - output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-msnsets_donor_1_2_3_preprocessed.rds
# ==============================================================================

library(Biostrings)
library(MSnSet.utils)
library(patchwork)
library(tidyverse)


# Import ------------------------------------------------------------------

source("data/single_islets_project_functions.R")

(load("output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData"))

# Load UniProt database
up_db <- read_tsv("data/uniprot_database_observed_proteins.tsv")

# Load FASTA for theoretical peptide calculation
fasta <- readAAStringSet(
	"data/Results_7_combined/2024-10-09-decoys-contam-uniprotkb_Human_AND_reviewed_true_AND_m_2024_10_09.fasta.fas"
)


# Process -----------------------------------------------------------------

# Calculate theoretical peptide count per protein (tryptic peptides 7-50 aa)
n_peptides <- as.character(fasta) %>%
	enframe() %>%
	filter(!grepl("^rev", name)) %>%
	mutate(UniProtAcc = str_extract(name, "(?<=\\|).+(?=\\|)"), seq = value) %>%
	select(-name) %>%
	rowid_to_column() %>%
	separate_longer_delim(cols = value, delim = stringr::regex("(?<=K|R)")) %>%
	mutate(pep_len = str_length(value)) %>%
	filter(pep_len >= 7, pep_len <= 50) %>%
	mutate(nPep = n(), .by = "rowid") %>%
	pull(nPep, UniProtAcc)


# Create output directory
dir.create(
	"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ",
	recursive = T,
	showWarnings = F
)

# Select Stage 1 T1D donors only (donor 3 in position 6)
msnsets <- msnsets[c(1, 2, 6)]

donor <- 1

for (i in names(msnsets)) {
	# Aggregate to protein level (sum intensities)
	msnsets[[i]] <- rrollup(msnsets[[i]], "ProteinName", "/", verbose = F, algorithm = "sum")

	# Calculate iBAQ: divide by theoretical peptide count
	uniprot_acc <- str_extract(featureNames(msnsets[[i]]), "(?<=\\|).+(?=\\|)")
	n_peptides_lp <- n_peptides[uniprot_acc]
	exprs(msnsets[[i]]) <- sweep(exprs(msnsets[[i]]), 1, STATS = n_peptides_lp, FUN = "/")

	# Normalize by column sum (fractional abundance)
	col_sum <- colSums(exprs(msnsets[[i]]), na.rm = T)
	exprs(msnsets[[i]]) <- sweep(exprs(msnsets[[i]]), 2, STATS = col_sum, FUN = "/")

	# Log2 transform
	exprs(msnsets[[i]]) <- log2(exprs(msnsets[[i]]))
	exprs(msnsets[[i]])[!is.finite(exprs(msnsets[[i]]))] <- NA

	# Parse UniProt identifier components
	fData(msnsets[[i]]) <- fData(msnsets[[i]]) %>%
		separate(col = ProteinName, into = c(NA, "UniProtAcc", "uniprot_id"), sep = "\\|")

	# Retain human proteins only
	msnsets[[i]] <- msnsets[[i]][grepl("_HUMAN", fData(msnsets[[i]])$uniprot_id), ]

	# Add gene name annotations
	fData(msnsets[[i]]) <- fData(msnsets[[i]]) %>%
		rownames_to_column("rowname") %>%
		left_join(up_db, by = c("UniProtAcc" = "From"), relationship = "one-to-one") %>%
		rename(gene_name = `Gene.Names..primary.`) %>%
		select(-Entry) %>%
		column_to_rownames("rowname")

	# Generate NA distribution plots
	png(
		filename = glue::glue(
			"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-Donor_{donor}_{i}_naplot.png"
		),
		width = 8,
		height = 8,
		units = "in",
		res = 200
	)
	naplot(msnsets[[i]])
	dev.off()

	# Identify outlier samples by missingness threshold
	mis_thresh <- 0.2
	sample_missingness <- apply(exprs(msnsets[[i]]), 2, \(x) {
		sum(is.na(x)) / length(x)
	})
	msnsets[[i]]$Outlier <- factor(sample_missingness > mis_thresh, levels = c("TRUE", "FALSE"))
	p5 <- plot_pca(msnsets[[i]], phenotype = "Outlier", show_ellipse = F)
	ggsave(
		filename = glue::glue(
			"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-Donor_{donor}_{i}_PCA_sample_outliers.png"
		),
		plot = p5,
		width = 4,
		height = 3,
		dpi = 600
	)

	# Exclude outlier samples
	msnsets[[i]] <- msnsets[[i]][, msnsets[[i]]$Outlier == FALSE]
	msnsets[[i]]$Outlier <- NULL

	# Pre-correction intensity distribution
	df_plot1 <- as.data.frame.MSnSet(msnsets[[i]]) %>%
		select(-make.names(setdiff(varLabels(msnsets[[i]]), "Block"))) %>%
		rownames_to_column("sample") %>%
		pivot_longer(cols = !matches("Block|sample"), names_to = "name", values_to = "value") %>%
		arrange(Block) %>%
		group_by(Block) %>%
		mutate(mean = mean(value, na.rm = T))

	p1 <- ggplot(df_plot1) +
		geom_boxplot(aes(x = fct_inorder(sample), y = value, color = Block)) +
		theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
		ggtitle(glue::glue("Donor {i} - Before Batch Correction")) +
		labs(color = "Pancreas Subregion") +
		xlab("Sample") +
		ylab("Intensity")

	ggsave(
		filename = glue::glue(
			"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-Donor_{donor}_{i}_before_batch_correction.png"
		),
		plot = p1,
		device = "png",
		scale = 1
	)

	# Pre-correction PCA visualization
	m_plot <- rename_block_for_plot(i)
	p3 <- plot_pca(m_plot, phenotype = "Block")
	ggsave(
		filename = glue::glue(
			"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-Donor_{donor}_{i}_before_batch_correction_PCA.png"
		),
		plot = p3,
		width = 4,
		height = 3,
		dpi = 600
	)

	# Apply tissue block batch correction
	msnsets[[i]] <- correct_batch_effect_NA(msnsets[[i]], "Block")

	# Post-correction intensity distribution
	df_plot2 <- as.data.frame.MSnSet(msnsets[[i]]) %>%
		select(-make.names(setdiff(varLabels(msnsets[[i]]), "Block"))) %>%
		rownames_to_column("sample") %>%
		pivot_longer(cols = !matches("Block|sample"), names_to = "name", values_to = "value") %>%
		arrange(Block) %>%
		group_by(Block) %>%
		mutate(mean = mean(value, na.rm = T))

	p2 <- ggplot(df_plot2) +
		geom_boxplot(aes(x = fct_inorder(sample), y = value, color = Block)) +
		theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
		ggtitle(glue::glue("Donor {i} - After Batch Correction")) +
		labs(color = "Pancreas Subregion") +
		xlab("Sample") +
		ylab("Intensity")

	ggsave(
		filename = glue::glue(
			"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-Donor_{donor}_{i}_after_batch_correction.png"
		),
		plot = p2,
		device = "png",
		scale = 1
	)

	# Post-correction PCA visualization
	m_plot2 <- rename_block_for_plot(i)
	p4 <- plot_pca(m_plot2, phenotype = "Block")
	ggsave(
		filename = glue::glue(
			"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-Donor_{donor}_{i}_after_batch_correction_PCA.png"
		),
		plot = p4,
		width = 4,
		height = 3,
		dpi = 600
	)

	donor <- donor + 1
}


saveRDS(
	msnsets,
	file = "output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-msnsets_donor_1_2_3_preprocessed.rds"
)
