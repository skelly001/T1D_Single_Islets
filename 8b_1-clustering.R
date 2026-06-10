# ==============================================================================
# WGCNA Co-Expression Network Clustering Analysis
# ==============================================================================
# Script: 8b_1-clustering.R
# Description: Performs refined clustering analysis on mAAb+ donors using
#              a Topological Overlap Matrix (TOM) of type signed Nowick 2 and
#              hierarchical clustering. Tests multiple soft-thresholding powers
#              and cluster numbers, validates clustering quality using GO enrichment
#              via clustGO function, and selects optimal parameters based on GO scores.
#
# Input: - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
# Output: - output/RD8-clustering/RD8b_1-clustering_results.RData
#         - output/RD8-clustering/RD8b_1-cluster_membership.xlsx
# ==============================================================================

library(MSnbase)
library(WGCNA)
library(tidyverse)
library(gprofiler2)


# Cluster comparison function --------------------------------------------

WGCNA_compare <- function(
	data,
	clusterNumbers = c(25, 50, 75, 100),
	powerNumber = c(9, 10, 11, 12, 13, 14),
	nameAlgorithm = c("TOM_signed_Nowick2")
) {
	dfList <- list()
	print('Initializing...')
	data <- t(scale(t(data), center = F))
	for (algorithm in nameAlgorithm) {
		for (p in powerNumber) {
			TOM <- TOMsimilarityFromExpr(
				t(data),
				weights = NULL,
				corType = "bicor",
				networkType = "signed",
				power = p,
				TOMType = "signed Nowick 2",
				TOMDenom = "min",
				maxPOutliers = 1,
				quickCor = 0,
				pearsonFallback = "individual",
				cosineCorrelation = F,
				replaceMissingAdjacencies = FALSE,
				suppressTOMForZeroAdjacencies = FALSE,
				suppressNegativeTOM = T,
				useInternalMatrixAlgebra = FALSE,
				nThreads = 0,
				verbose = 1,
				indent = 0
			)
			rownames(TOM) <- rownames(data)
			colnames(TOM) <- rownames(data)
			TOM_cut <- hclust(as.dist(1 - TOM), method = 'ward.D')
			for (clusters in clusterNumbers) {
				if (algorithm == 'TOM_signed_Nowick2') {
					print('Running TOM_signed_Nowick_average')

					result <- dendextend::cutree(TOM_cut, k = clusters, order_clusters_as_data = F)

					df <- data.frame(
						Gene = names(result),
						Cluster = result
					)

					# Assign a unique Identifier
					df_name <- paste(algorithm, "_Power", p, "_wardD_clusters_", clusters, sep = "")

					# Append to Results List
					dfList[[df_name]] <- df
				}
			}
		}
	}
	dfList <- dfList[sapply(dfList, nrow) > 0]
	return(dfList)
}


# Fast Correlation Function -----------------------------------------------

fast_cor <- function(data, use = "pairwise.complete.obs", method = "spearman") {
	{
		n = NULL
		n <- psych::pairwiseCount(data) - 2 # Create matrix of degrees freedom
	}
	if (method == "spearman") {
		data <- base::apply(X = data, MARGIN = 2, data.table::frankv)
	}
	r = NULL
	p1 = NULL
	p2 = NULL
	r <- coop::pcor(x = data, use = use)
	{
		t <- sqrt(n) * r / sqrt(1 - r^2) # T-test
		p1 <- stats::pt(t, n) # One side p-value
		p2 <- stats::pt(t, n, lower.tail = FALSE)
	}
	flt.Corr.Matrix <- function(cormat, pmat1 = NULL, pmat2 = NULL, df = NULL) {
		ut <- base::upper.tri(cormat)
		flt_data <- data.frame(
			Var1 = base::rownames(cormat)[base::row(cormat)[ut]],
			Var2 = base::rownames(cormat)[base::col(cormat)[ut]],
			cor = cormat[ut]
		)
		if (!is.null(pmat1)) {
			flt_data$p1 <- pmat1[ut]
		}
		if (!is.null(pmat2)) {
			flt_data$p2 <- pmat2[ut]
		}
		if (!is.null(df)) {
			flt_data$df <- df[ut]
		}
		return(flt_data)
	}
	{
		result <- flt.Corr.Matrix(cormat = r, df = n, pmat1 = p1, pmat2 = p2)
		result <- base::transform(result, p = base::pmin(p1, p2) * 2)
		result$FDR <- stats::p.adjust(result$p, method = "BH")
		result <- base::subset(result, select = -c(p1, p2))
		result <- list(result, r)
	}
	return(result)
}


# Import -----------------------------------------------------------------

source("data/single_islets_project_functions.R")
source("8a_1-clustGO_Function_GOScore.R")

# Load protein abundance data with centroids
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")


# Process ----------------------------------------------------------------

dir.create("output/RD8-clustering")

# Filter, impute, and subset to mAAb+ donors
m <- completeness_filt(m, 0.5)
m <- impute_msnset(m)
m <- m[, pData(m)$class == "case"]

# Define GO background as all genes in filtered dataset
bg <- rownames(exprs(m))

## Run WGCNA Clustering ----
WGCNA_Clusters <- WGCNA_compare(
	exprs(m),
	clusterNumbers = c(75),
	nameAlgorithm = c("TOM_signed_Nowick2"),
	powerNumber = c(12)
)

WGCNA_Summary <- clustGO(
	WGCNA_Clusters,
	organism = "hsapiens",
	Bfactor = 0.5,
	S_type = "Summary",
	bg = bg
)

# Compute Pearson correlation matrix
correlations_combined <- fast_cor(t(exprs(m)), method = "pearson")
correlations_df <- correlations_combined[[1]]
correlations_matrix <- correlations_combined[[2]]

# Save
save(
	WGCNA_Clusters,
	WGCNA_Summary,
	correlations_combined,
	file = "output/RD8-clustering/RD8b_1-clustering_results.RData"
)

# Export cluster membership
openxlsx::write.xlsx(WGCNA_Summary[[4]], "output/RD8-clustering/RD8b_1-cluster_membership.xlsx")
