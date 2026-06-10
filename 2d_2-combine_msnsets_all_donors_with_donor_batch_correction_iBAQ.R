# ==============================================================================
# Combine All Donor iBAQ MSnSets With Donor Batch Correction
# ==============================================================================
# Script: 2d_2-combine_msnsets_all_donors_with_donor_batch_correction_iBAQ.R
# Description: Combines preprocessed iBAQ quantification data from all six donors
#              (mAAb+ and control) into a single MSnSet object with donor-level
#              batch correction. Adds case/control status and class labels based
#              on insulin (INS) and CD3 status.
#
# Input: - output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-msnsets_donor_1_2_3_preprocessed.rds
#        - output/RD2-preprocessing/RD2c_2-control_donors_preprocessing_iBAQ/RD2c_2-msnsets_donor_4_5_6_preprocessed.rds
# Output: - output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_2-msnset_all_donors_BC_iBAQ.rds
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

all <- vctrs::vec_c(case, control)


# Merge MSnSets into single object
all2 <- purrr::reduce(all, MSnbase::combine)


# Apply donor-level batch correction
all3 <- correct_batch_effect_NA(all2, batch_name = "donor")


# Assign feature names as gene_UniProt format
featureNames(all3) <- paste(fData(all3)$gene_name, fData(all3)$UniProtAcc, sep = "...")


# Set data types and add pData
pData(all3)$is.case <- pData(all3)$donor %in% c("6450", "6521", "6267")
pData(all3)$donor <- factor(
	pData(all3)$donor,
	levels = c("6450", "6521", "6267", "6178", "6440", "6539")
)
pData(all3)$class <- ifelse(pData(all3)$is.case == TRUE, "case", "control")
pData(all3)$class <- factor(pData(all3)$class, levels = c("control", "case"))
pData(all3)$class_ins <- paste(pData(all3)$class, pData(all3)$INS, sep = "_")
pData(all3)$INS <- factor(pData(all3)$INS, levels = c("Negative", "Positive"))
pData(all3)$CD3 <- factor(pData(all3)$CD3, levels = c("Negative", "Positive"))
fData(all3)$feature <- rownames(fData(all3))
pData(all3)$MS_block <- case_when(
	pData(all3)$donor == "6450" ~ "MSB1",
	pData(all3)$donor == "6521" ~ "MSB2",
	.default = "MSB3"
) %>%
	factor(., levels = c("MSB1", "MSB2", "MSB3"))


# Remove NaN values, all-NA features, and duplicate gene names
exprs(all3)[is.nan(exprs(all3))] <- NA
all3 <- all3[
	apply(exprs(all3), 1, \(x) {
		!all(is.na(x))
	}),
]
all3 <- all3[!is.na(fData(all3)$gene_name), ]
fData(all3) <- fData(all3) %>%
	rownames_to_column() %>%
	group_by(gene_name) %>%
	mutate(count = n()) %>%
	column_to_rownames()
all3 <- all3[fData(all3)$count == 1, ]
featureNames(all3) <- fData(all3)$gene_name
fData(all3)$count <- NULL
fData(all3)$feature <- NULL


saveRDS(
	all3,
	"output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_2-msnset_all_donors_BC_iBAQ.rds"
)
