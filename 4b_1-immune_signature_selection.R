# ==============================================================================
# Immune Signature Feature Selection Using Machine Learning
# ==============================================================================
# Script: 4b_1-immune_signature_selection.R
# Description: Refines the Islet Immune Response Signature (IIRS) using machine
#              learning-based feature selection. Employs random forest models with
#              recursive feature elimination to identify the most discriminative
#              proteins distinguishing INS+/CD3+ islets from other classes.
#              Aggregates results across multiple runs using Robust Rank
#              Aggregation.
#
# Input: - output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds
#        - output/RD4-islet_immune_response_signature/RD4a_1-WGCNA_module_immune_signature_candidate_proteins.xlsx
# Output: - output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds
#         - output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)
library(mlr3verse)
library(mlr3fselect)
library(mlr3tuning)
library(mlr3viz)
library(future)
library(RobustRankAggreg)
plan(multisession, workers = 6)


# Import ------------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load batch-corrected protein abundance data
m <- readRDS("output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds")

# Load WGCNA-derived immune signature candidates
im_prot <- openxlsx::read.xlsx(
	"output/RD4-islet_immune_response_signature/RD4a_1-WGCNA_module_immune_signature_candidate_proteins.xlsx"
)

# Process -----------------------------------------------------------------

# Subset to immune module candidate proteins
m2 <- m[fData(m)$UniProtAcc %in% im_prot$UniProtAcc, ]

# Require 50% completeness per donor
m3 <- completeness_filt(m2, 0.5)


# Set feature names as gene_UniProt format
featureNames(m3) <- paste(fData(m3)$gene_name, fData(m3)$UniProtAcc, sep = "...")


# Feature Select ----------------------------------------------------------

data <- as.data.frame.MSnSet(m3) %>%
	select(-base::setdiff(varLabels(m3), c("class"))) %>%
	rownames_to_column("ID")


if (!file.exists("output/RD4-islet_immune_response_signature/RD4b_1-final_ML_models.rds")) {
	set.seed(0)

	# Task
	task = TaskClassif$new("data_task", backend = data, target = "class")

	# Learner
	learner = lrn(
		"classif.ranger",
		predict_type = "prob",
		importance = "permutation",
		num.threads = 6,
		num.trees = 1000
	)

	# Update column roles
	task$set_col_roles(
		cols = "ID",
		add_to = "name",
		remove_from = "feature"
	)
	task$set_col_roles(
		cols = "class",
		add_to = "stratum"
	)

	# Inner and outer resampling
	inner = rsmp("cv", folds = 5)
	outer = rsmp("cv", folds = 10)

	# Feature selector
	fselector = fs("rfe", n_features = 60, feature_fraction = 0.9)

	# Hyperparameter tuner
	tuner = tnr("random_search", batch_size = 5)
	search_space = ps(
		mtry = p_int(1, 20),
		min.node.size = p_int(5, 20)
	)

	# Combine tuning into learner
	auto_tuned = AutoTuner$new(
		learner = learner,
		resampling = inner,
		measure = msr("classif.auc"),
		tuner = tuner,
		search_space = search_space,
		terminator = trm("evals", n_evals = 10),
		store_tuning_instance = TRUE
	)

	# Nested feature selection with tuning
	rr = fselect_nested(
		fselector = fselector,
		task = task,
		learner = auto_tuned,
		inner_resampling = inner,
		outer_resampling = outer,
		measure = msr("classif.auc"),
		terminator = trm("none"),
		store_fselect_instance = TRUE
	)

	saveRDS(rr, "output/RD4-islet_immune_response_signature/RD4b_1-final_ML_models.rds")
} else {
	rr <- readRDS("output/RD4-islet_immune_response_signature/RD4b_1-final_ML_models.rds")
}


# Extract selected features from each outer CV fold
afs_list = rr$learners
rr_feat <- map(afs_list, \(x) {
	{
		as.data.frame(x$archive$data) %>%
			slice_tail(n = 1) %>%
			pull(importance) %>%
			.[[1]]
	} %>%
		names()
})


# Aggregate feature rankings using Robust Rank Aggregation
ag_rank <- aggregateRanks(glist = rr_feat, method = "RRA", exact = T)


# Select top 40 proteins by aggregated rank score
fin_feat <- ag_rank %>%
	arrange(Score) %>%
	slice_head(n = 40) %>%
	pull(Name)

# Final IIRS (40 proteins)
im_sig_out <- tibble(
	rank = 1:40,
	gene_name = str_extract(fin_feat, ".+(?=\\.\\.\\.)"),
	UniProtAcc = str_extract(fin_feat, "(?<=\\.\\.\\.).+")
)

saveRDS(
	im_sig_out,
	"output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds"
)


# Export all ranked candidate proteins
im_sig_all_candidates <- ag_rank %>%
	arrange(Score) %>%
	mutate(
		rank = 1:nrow(.),
		gene_name = str_extract(Name, ".+(?=\\.\\.\\.)"),
		UniProtAcc = str_extract(Name, "(?<=\\.\\.\\.).+")
	) %>%
	select(rank, gene_name, UniProtAcc, P.Value = Score, -Name) %>%
	`rownames<-`(NULL)

openxlsx::write.xlsx(
	im_sig_all_candidates,
	"output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_feature_candidates.xlsx"
)


# Calculate IIRS centroid (median of signature proteins) per sample
m <- readRDS("output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds")
m_im <- m[fData(m)$UniProtAcc %in% im_sig_out$UniProtAcc, ]
pData(m)$im_sig_centroid <- apply(exprs(m_im), 2, \(x) {
	median(x, na.rm = TRUE)
})
saveRDS(m, "output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds")
