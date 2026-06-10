# ==============================================================================
# Beta Cell Profile Differential Expression Analysis
# ==============================================================================
# Script: 7a_1-beta_cell_profile_differential_analysis.R
# Description: Performs differential protein expression analysis for Beta Cell
#              Profile using limma. Applies 50% completeness filtering, standardizes
#              BCP centroid values (Z-score) within donors for CAMERA-PR, and performs
#              donor-specific linear modeling to identify proteins associated with
#              beta cell function.
#
# Input: - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
# Output: - output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv
# ==============================================================================

library(ComplexHeatmap)
library(MSnSet.utils)
library(tidyverse)


# Import -----------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load protein abundance data with centroids
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")


# Process ----------------------------------------------------------------

dir.create("output/RD7-beta_cell_profile_analysis")


# Require 50% completeness per donor
m2 <- completeness_filt(m, 0.5)

# Z-score BCP centroid within donors for comparable effect sizes
pData(m2) <- pData(m2) %>%
	mutate(.by = donor, im_zscore = scale(bcp_centroid)[, 1])


## Differential analysis ----

### Linear modeling across all islets ----

# Limma by donor
res <- lapply(levels(m2$donor), function(donor_i) {
	limma_gen(m2[, m2$donor == donor_i], "~ im_zscore", "im_zscore") %>%
		rownames_to_column("protein")
}) %>%
	setNames(levels(m2$donor)) %>%
	bind_rows(.id = "donor") %>%
	mutate(donor = factor(donor, levels = unique(donor))) %>%
	select(-B, -starts_with("SE", ignore.case = FALSE)) %>%
	mutate(ZScore = qnorm(P.Value / 2, lower.tail = FALSE) * sign(logFC))

table(res$donor, res$adj.P.Val < 0.05)

# Save
out <- res %>%
	pivot_wider(
		id_cols = protein,
		names_from = donor,
		values_from = c(logFC, ZScore, P.Value, adj.P.Val)
	)

write.csv(
	out,
	file = "output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv",
	quote = FALSE,
	row.names = FALSE
)


#############################################################################
# Reviewer #1 comment - is ECM causal or caused? (I dont think this can really be
# answered with our data. WQ wants to just look at correlation)

### ECM ttest ----

ECM_prot <- c("FAP", "DPP4", "COL4A1", "COL4A2", "COL4A3")
m2_ecm <- m2[fData(m2)$gene_name %in% ECM_prot, ]
m2_ecm$INS_low <- exprs(m2)["INS", ] < -2


### K-means group comparison ----
df <- as.data.frame.MSnSet(m2_ecm)

# K-means clustering within each donor
set.seed(0)
df <- df %>%
	group_by(donor) %>%
	mutate(
		bcp_cluster = {
			km_i <- kmeans(bcp_centroid, centers = 4, iter.max = 100, nstart = 100)
			# Relabel clusters by ascending center value
			rank_map <- rank(km_i$centers[, 1])
			rank_map[km_i$cluster]
		}
	) %>%
	ungroup()

ggplot(df, aes(ordered(bcp_cluster), bcp_centroid)) +
	geom_boxplot() +
	facet_wrap(~donor)

# BCP-low (rank 1) vs low-intermediate (rank 2)
df_km <- df %>%
	filter(bcp_cluster %in% c(1, 2)) %>%
	mutate(group = if_else(bcp_cluster == 1, "BCP-low", "Low-intermediate"))
# %>%
# filter(!INS_low)

# Test all ECM proteins between k-means groups, within each donor
km_tests <- expand.grid(
	protein = ECM_prot,
	donor = unique(df_km$donor),
	stringsAsFactors = FALSE
) %>%
	pmap_dfr(function(protein, donor) {
		d <- df_km %>% filter(donor == !!donor)
		test <- t.test(reformulate("group", protein), data = d)
		tibble(protein = protein, donor = donor, W = test$statistic, p_value = test$p.value)
	})

print(km_tests, n = 1000)

# Boxplots for all ECM proteins (all 4 clusters)
df %>%
	pivot_longer(cols = all_of(ECM_prot), names_to = "protein", values_to = "intensity") %>%
	ggplot(aes(x = ordered(bcp_cluster), y = intensity)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter(width = 0.2, alpha = 0.2, size = 1.5) +
	facet_grid(donor ~ protein, scales = "free_y") +
	labs(x = "K-means cluster", y = "Intensity (log2)") +
	theme_grey(16)

# Boxplots for all ECM proteins (BCP-low vs Low-intermediate)
df_km %>%
	pivot_longer(cols = all_of(ECM_prot), names_to = "protein", values_to = "intensity") %>%
	ggplot(aes(x = group, y = intensity)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter(width = 0.2, alpha = 0.2, size = 1.5) +
	facet_wrap(~protein, scales = "free_y") +
	labs(x = NULL, y = "Intensity (log2)") +
	theme_grey(16)


### Quantile-based group comparison ----

df <- df %>%
	group_by(donor) %>%
	mutate(
		# .by = donor,
		bcp_group = cut(
			bcp_centroid,
			breaks = quantile(bcp_centroid, probs = c(0, 0.25, 0.5, 0.75, 1)),
			labels = c("Q1-low", "Q2", "Q3", "Q4-high"),
			include.lowest = TRUE
		)
	)

table(df$donor, df$bcp_group)

df_q <- df %>%
	filter(bcp_group %in% c("Q1-low", "Q2"))

# Q1/Q2 may be skewed (so wilcox)
df_q %>% filter(bcp_group == "Q1-low") %>% pull(DPP4) %>% hist()
df_q %>% filter(bcp_group == "Q2") %>% pull(DPP4) %>% hist()

# Test all ECM proteins between quantile groups, within each donor
q_tests <- expand.grid(protein = ECM_prot, donor = unique(df_q$donor), stringsAsFactors = FALSE) %>%
	pmap_dfr(function(protein, donor) {
		d <- df_q %>% filter(donor == !!donor)
		test <- t.test(reformulate("bcp_group", protein), data = d)
		tibble(protein = protein, donor = donor, t = test$statistic, p_value = test$p.value)
	})

print(q_tests, n = 1000)

# Boxplots for all ECM proteins (all clusters, T1D donors only)
t1d_donors <- c("6450", "6521", "6267")

# Create significance labels
q_sig <- q_tests %>%
	mutate(sig = if_else(p_value < 0.05, "p < 0.05", "n.s."))

df %>%
	filter(donor %in% t1d_donors) %>%
	pivot_longer(cols = all_of(ECM_prot), names_to = "protein", values_to = "intensity") %>%
	left_join(q_sig, by = c("donor", "protein")) %>%
	mutate(
		fill_group = if_else(
			bcp_group %in% c("Q1-low", "Q2") & sig == "p < 0.05",
			"p < 0.05",
			"n.s."
		)
	) %>%
	ggplot(aes(x = bcp_group, y = intensity, fill = fill_group)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter(width = 0.2, alpha = 0.4) +
	scale_fill_manual(values = c("p < 0.05" = "#F8766D", "n.s." = "grey80")) +
	facet_grid(donor ~ protein, scales = "free_y") +
	labs(x = "BCP", y = "Intensity (log2)", fill = NULL) +
	theme_grey(16)

# Boxplots for all ECM proteins (Q1-low vs Q2, T1D donors only)
df_q %>%
	filter(donor %in% t1d_donors) %>%
	pivot_longer(cols = all_of(ECM_prot), names_to = "protein", values_to = "intensity") %>%
	left_join(q_sig, by = c("donor", "protein")) %>%
	ggplot(aes(x = bcp_group, y = intensity, fill = sig)) +
	geom_boxplot(outlier.shape = NA) +
	geom_jitter(width = 0.2, alpha = 0.4) +
	scale_fill_manual(values = c("p < 0.05" = "#F8766D", "n.s." = "grey80"), na.value = "grey80") +
	facet_grid(donor ~ protein, scales = "free_y") +
	labs(x = "BCP", y = "Intensity (log2)", fill = NULL) +
	theme_grey(16)

#### INS low

### BCP vs FAP ----

ep_data <- as.data.frame.MSnSet(m2) %>%
	mutate(INS_low = INS < -2) %>%
	mutate(bcp_centroid = scale(bcp_centroid), .by = donor)

ggplot(ep_data, aes(bcp_centroid, DPP4)) +
	facet_wrap("donor") +
	geom_point() +
	geom_smooth(method = "lm") +
	ggpubr::stat_cor() +
	theme_grey(base_size = 16)

ggplot(ep_data, aes(bcp_centroid, FAP, color = INS_low)) +
	geom_boxplot() +
	ggpubr::stat_cor() +
	theme_grey(base_size = 16)
