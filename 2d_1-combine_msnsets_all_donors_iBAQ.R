# ==============================================================================
# Combine All Donor iBAQ MSnSets Without Batch Correction
# ==============================================================================
# Script: 2d_1-combine_msnsets_all_donors_iBAQ.R
# Description: Combines preprocessed iBAQ quantification data from all six donors
#              (Stage 1 T1D and control) into a single MSnSet object. Adds case/control
#              status and class labels based on insulin (INS) and CD3 status without
#              donor-level batch correction.
#
# Input: - output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-msnsets_donor_1_2_3_preprocessed.rds
#        - output/RD2-preprocessing/RD2c_2-control_donors_preprocessing_iBAQ/RD2c_2-msnsets_donor_4_5_6_preprocessed.rds
# Output: - output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_1-msnset_all_donors_iBAQ.rds
# ==============================================================================

library(org.Hs.eg.db)
library(MSnSet.utils)
library(tidyverse)


# Import ------------------------------------------------------------------

case <- readRDS(
	"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-msnsets_donor_1_2_3_preprocessed.rds"
)

control <- readRDS(
	"output/RD2-preprocessing/RD2c_2-control_donors_preprocessing_iBAQ/RD2c_2-msnsets_donor_4_5_6_preprocessed.rds"
)


# Process -----------------------------------------------------------------

dir.create("output/RD2-preprocessing/RD2d-prepared_data_iBAQ", showWarnings = F)


all <- vctrs::vec_c(case, control)


# Merge MSnSets into single object
all2 <- purrr::reduce(all, MSnbase::combine)

# Assign feature names as gene_UniProt format
featureNames(all2) <- paste(fData(all2)$gene_name, fData(all2)$UniProtAcc, sep = "...")


# Set data types and add pData
pData(all2)$is.case <- pData(all2)$donor %in% c("6450", "6521", "6267")
pData(all2)$donor <- factor(
	pData(all2)$donor,
	levels = c("6450", "6521", "6267", "6178", "6440", "6539")
)
pData(all2)$class <- ifelse(pData(all2)$is.case == TRUE, "case", "control")
pData(all2)$class <- factor(pData(all2)$class, levels = c("control", "case"))
pData(all2)$class_ins <- paste(pData(all2)$class, pData(all2)$INS, sep = "_")
pData(all2)$INS <- factor(pData(all2)$INS, levels = c("Negative", "Positive"))
pData(all2)$CD3 <- factor(pData(all2)$CD3, levels = c("Negative", "Positive"))
fData(all2)$feature <- rownames(fData(all2))
pData(all2)$MS_block <- case_when(
	pData(all2)$donor == "6450" ~ "MSB1",
	pData(all2)$donor == "6521" ~ "MSB2",
	.default = "MSB3"
) %>%
	factor(., levels = c("MSB1", "MSB2", "MSB3"))


# Remove NaN values, all-NA features, and duplicate gene names
exprs(all2)[is.nan(exprs(all2))] <- NA
all2 <- all2[
	apply(exprs(all2), 1, \(x) {
		!all(is.na(x))
	}),
]
all2 <- all2[!is.na(fData(all2)$gene_name), ]
fData(all2) <- fData(all2) %>%
	rownames_to_column() %>%
	group_by(gene_name) %>%
	mutate(count = n()) %>%
	column_to_rownames()
all2 <- all2[fData(all2)$count == 1, ]
featureNames(all2) <- fData(all2)$gene_name
fData(all2)$count <- NULL
fData(all2)$feature <- NULL


saveRDS(all2, "output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_1-msnset_all_donors_iBAQ.rds")
