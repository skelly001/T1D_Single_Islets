# ==============================================================================
# Immune Signature Candidate Protein Selection
# ==============================================================================
# Script: 4a_1-immune_signature_candidate_proteins.R
# Description: Identifies candidate proteins for the Islet Immune Response
#              Signature (IIRS) by extracting proteins from immune-related WGCNA
#              modules across the three mAAb+ donors. Consolidates module
#              membership results and exports candidate proteins for further
#              refinement.
#
# Input: - WGCNA module membership files from output/RD3-WGCNA/ for donors 6450, 6521, 6267
#        - output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds
# Output: - output/RD4-islet_immune_response_signature/RD4a_1-WGCNA_module_immune_signature_candidate_proteins.xlsx
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)


# Import ------------------------------------------------------------------

m <- readRDS("output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds")

# Load WGCNA module membership files for case donors
file_list_case <- list.files(
	"output/RD3-WGCNA",
	pattern = "WGCNA_module_membership",
	full.names = T,
	recursive = T
)

file_list_case <- tibble(paths = file_list_case) %>%
	mutate(path_order = str_extract(paths, "(?<=/RD3)[[:alpha:]]")) %>%
	arrange(path_order) %>%
	pull(paths)

res_case <- map(file_list_case, readxl::read_xlsx)


# Process -----------------------------------------------------------------

dir.create("output/RD4-islet_immune_response_signature")

# Extract donor identifiers from file paths
file_names_case <- str_extract(file_list_case, "donor_\\d{4}")
names(res_case) <- file_names_case

# Define immune-related modules per donor (manually curated from WGCNA results)
sig_clust <- list(
	donor_6450 = c("royalblue", "lightyellow"),
	donor_6521 = "green",
	donor_6267 = "greenyellow"
)

# Filter proteins belonging to immune-related modules
res2 <- map2(res_case, sig_clust, \(x, y) {
	filter(x, moduleColors %in% y) %>%
		dplyr::select(UniProtAcc:moduleColors) %>%
		distinct()
})

res3 <- res2 %>%
	enframe() %>%
	unnest("value") %>%
	rename(Donor = name) %>%
	mutate(modules = paste(Donor, moduleColors, sep = "-"))

openxlsx::write.xlsx(
	res3,
	file = "output/RD4-islet_immune_response_signature/RD4a_1-WGCNA_module_immune_signature_candidate_proteins.xlsx"
)
