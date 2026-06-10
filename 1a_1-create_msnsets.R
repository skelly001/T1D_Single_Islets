# ==============================================================================
# MSnSet Creation from MSstats Data
# ==============================================================================
# Script: 1a_1-create_msnsets.R
# Description: Creates MSnSet objects from MSstats quantification results for
#              all six donors. Processes protein abundance data, creates
#              wide-format matrices, integrates phenotype metadata, and generates
#              separate MSnSet objects for experimental samples and reference
#              pools.
#
# Input: - MSstats CSV files in data/Results_7_combined/
#        - data/sample_metadata.xlsx
# Output: - output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData
#         - output/RD1-raw_msnsets/RD1a_1-msnsets_ref_raw.RData
# ==============================================================================

library(MSnSet.utils) # remotes::install_github("PNNL-Comp-Mass-Spec/MSnSet.utils")
library(readxl)
library(tidyverse)
library(data.table)


# Abundance Data ----------------------------------------------------------

dir.create("output/RD1-raw_msnsets", recursive = TRUE, showWarnings = F)

# Donor mapping (Data Package ID -> Patient ID):
# MSnSet 1: pid = 6450, DP3802 (mAAb+)
# MSnSet 2: pid = 6521, DP4030 (mAAb+)
# MSnSet 3: pid = 6178, DP5312 (ND)
# MSnSet 4: pid = 6440, DP5313 (ND)
# MSnSet 5: pid = 6539, DP5314 (ND)
# MSnSet 6: pid = 6267, DP5315 (mAAb+)

pid <- c(
	"3802" = "6450",
	"4030" = "6521",
	"5312" = "6178",
	"5313" = "6440",
	"5314" = "6539",
	"5315" = "6267"
)

# Load MSstats quantification files
file_list <- list.files(
	path = "data/Results_7_combined",
	pattern = "MSstats.csv",
	full.names = T,
	recursive = T
)
msstats <- rbindlist(lapply(file_list, fread, showProgress = TRUE, data.table = TRUE))

# Extract DP identifier from condition string
msstats[, DP := str_extract(Condition, "(?<=^DP)\\d{4}")]

# Split data by DP
msstats <- base::split(msstats, by = "DP") %>%
	map(\(x) {
		select(x, -DP)
	})


# Generate sample names based on donor-specific naming conventions
for (i in names(msstats)) {
	if (i == "3802") {
		msstats[[i]][,
			MeasurementName := fcase(
				grepl("Reference", Condition) , paste0("Ref_", str_extract(Condition, "SW.*?_\\D\\d{1,2}$")) ,
				default = str_extract(Condition, "SW.*?_\\D\\d{1,2}$")
			)
		]
	} else if (i == "4030") {
		msstats[[i]][,
			MeasurementName := fcase(
				grepl("Reference", Condition) , paste0("Ref", str_extract(Condition, "(?<=Reference).*")) ,
				default = str_extract(Condition, "SW.*?_\\D\\d{1,2}$")
			)
		]
	} else if (i %in% c("5312", "5313", "5314")) {
		msstats[[i]][, MeasurementName := str_extract(Condition, "Chip.*\\D\\d{1,2}$")]
	} else if (i == "5315") {
		msstats[[i]][, MeasurementName := sub(".*(Chip.*)(\\D\\d{1,2})$", "\\1_\\2", Condition)]
	} else if (i == "5506") {
		msstats[[i]][, MeasurementName := paste0("Ref", str_extract(Condition, "(?<=Reference).*"))]
	}
}

# Convert to data frames
msstats <- map(msstats, as.data.frame)


# Separate donor 3-6 reference samples
msnsets_ref <- vector("list", 3)
names(msnsets_ref) <- c("6450", "6521", "6178_6440_6539_6267")
msnsets_ref[["6178_6440_6539_6267"]] <- msstats[["5506"]]
msstats[["5506"]] <- NULL

# Rename list elements with patient IDs
names(msstats) <- pid[names(msstats)]


# Create MSnSets ----------------------------------------------------------

# Sum function returning NA when all values are NA
sum_na_rm <- \(x) {
	if (all(is.na(x))) {
		NA
	} else {
		sum(x, na.rm = TRUE)
	}
}

wide_data <- \(x) {
	x2 <- x %>%
		filter(!grepl("^rev", ProteinName)) %>%
		mutate(FeatureName = paste(ProteinName, PeptideSequence, PrecursorCharge, sep = "-")) %>%
		dplyr::select(
			FeatureName,
			ProteinName,
			MeasurementName,
			PeptideSequence,
			PrecursorCharge,
			Intensity
		) %>%
		pivot_wider(
			names_from = MeasurementName,
			values_from = Intensity,
			values_fn = sum_na_rm
		) %>%
		as.data.frame()

	rownames(x2) <- x2$FeatureName

	return(x2)
}

# Pivot to wide format
msstats <- map(msstats, wide_data)


# Pivot donor 3-6 reference data to wide format
msnsets_ref[["6178_6440_6539_6267"]] <- wide_data(msnsets_ref[["6178_6440_6539_6267"]])


# Construct MSnSet objects
msnsets <- map(msstats, \(x) {
	e_data <- select(x, -c(FeatureName:PrecursorCharge)) %>% as.matrix()
	f_data <- select(x, c(FeatureName:PrecursorCharge))
	MSnSet(exprs = e_data, fData = f_data)
})
rm(msstats)


# Separate reference samples from donors 1-2
donors_to_split <- c("6450", "6521")
for (i in donors_to_split) {
	msnsets_ref[[i]] <- msnsets[[i]][, grepl("Ref", sampleNames(msnsets[[i]]))]
	msnsets[[i]] <- msnsets[[i]][, !grepl("Ref", sampleNames(msnsets[[i]]))]
}

# Construct MSnSet for donor 3-6 references
e_data <- select(msnsets_ref[["6178_6440_6539_6267"]], -c(FeatureName:PrecursorCharge)) %>%
	as.matrix()
f_data <- select(msnsets_ref[["6178_6440_6539_6267"]], c(FeatureName:PrecursorCharge))
msnsets_ref[["6178_6440_6539_6267"]] <- MSnSet(exprs = e_data, fData = f_data)


# Phenotype Data ----------------------------------------------------------

# Load sample metadata
pdata <- readxl::read_xlsx("data/sample_metadata.xlsx", sheet = "experiment_samples")

pdata <- pdata %>%
	rename(
		Block = tissue_block_id,
		file_name = `file_name\r\n(file name for processing; DP + raw_dataset_file_name + PID)`
	)


## Donor 1 (6450) ----

# Merge phenotype data with MSnSet
pData(msnsets[["6450"]]) <- pData(msnsets[["6450"]]) %>%
	rownames_to_column("dataset_id") %>%
	left_join(pdata, by = "dataset_id") %>%
	column_to_rownames("dataset_id")

# Exclude sample with high missingness
msnsets[["6450"]] <- msnsets[["6450"]][, !grepl("SW190842_C5", sampleNames(msnsets[["6450"]]))]


## Donor 2 (6521) ----

# Merge phenotype data with MSnSet
pData(msnsets[["6521"]]) <- pData(msnsets[["6521"]]) %>%
	rownames_to_column("dataset_id") %>%
	left_join(pdata, by = "dataset_id") %>%
	column_to_rownames("dataset_id")

# Exclude sample with high missingness
msnsets[["6521"]] <- msnsets[["6521"]][, !grepl("SW210407_A3", colnames(msnsets[["6521"]]))]


## Donors 3-5 (6178, 6440, 6539) ----

donors <- c("6178", "6440", "6539")
for (i in donors) {
	pdata_lp <- pData(msnsets[[i]]) %>%
		rownames_to_column("dataset_id") %>%
		left_join(pdata, by = "dataset_id") %>%
		column_to_rownames("dataset_id")

	sampleNames(msnsets[[i]]) <- rownames(pdata_lp)
	pData(msnsets[[i]]) <- pdata_lp
}


## Donor 6 (6267) ----

# Merge phenotype data with MSnSet
pdata_d6 <- pData(msnsets[["6267"]]) %>%
	rownames_to_column("dataset_id") %>%
	left_join(pdata, by = "dataset_id") %>%
	column_to_rownames("dataset_id")

sampleNames(msnsets[["6267"]]) <- rownames(pdata_d6)
pData(msnsets[["6267"]]) <- pdata_d6


# Standardize Sample Names ----
msnsets <- map(msnsets, \(x) {
	sampleNames(x) <- pData(x)$sample_name
	return(x)
})


# Save ----
save(msnsets, file = "output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData")
save(msnsets_ref, file = "output/RD1-raw_msnsets/RD1a_1-msnsets_ref_raw.RData")
