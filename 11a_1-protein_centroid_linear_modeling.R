# ==============================================================================
# Protein–Centroid Linear Modeling
# ==============================================================================
# Script: 11a_1-protein_centroid_linear_modeling.R
# Description: Models selected protein abundances against beta cell profile (BCP)
#              or islet immune response signature (IIRS) centroids. Loads protein
#              abundance data with centroids, removes islet selection outliers, and
#              generates per-donor scatterplots with fitted linear models, saving
#              one PNG per protein. Covers hyaluronan-related proteins (ITIH5, HEXA,
#              HEXB) vs IIRS and extracellular matrix proteins (FAP, DPP4, QSOX1,
#              FBLN7) vs BCP.
#
# Input: - output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds
# Output: - output/RD11-protein_centroid_linear_modeling/RD11a_1-ITIH5_vs_IIRS.png
#         - output/RD11-protein_centroid_linear_modeling/RD11a_1-HEXA_vs_IIRS.png
#         - output/RD11-protein_centroid_linear_modeling/RD11a_1-HEXB_vs_IIRS.png
#         - output/RD11-protein_centroid_linear_modeling/RD11a_1-FAP_vs_BCP.png
#         - output/RD11-protein_centroid_linear_modeling/RD11a_1-DPP4_vs_BCP.png
#         - output/RD11-protein_centroid_linear_modeling/RD11a_1-QSOX1_vs_BCP.png
#         - output/RD11-protein_centroid_linear_modeling/RD11a_1-FBLN7_vs_BCP.png
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)
library(patchwork)

# Import -----------------------------------------------------------------

# Protein abundance data with centroids
m <- readRDS("output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds")


# Process -----------------------------------------------------------------

# Output location
base_path <- "output/RD11-protein_centroid_linear_modeling"
fs::dir_create(base_path)


# Evaluate islet selection outliers
ep_data <- as.data.frame.MSnSet(m) %>%
	group_by(donor) %>%
	mutate(BCP_GCG_low = (bcp_centroid < -1) & (GCG < -2))


# Remove outliers
plot_data <- ep_data %>%
	ungroup() %>%
	filter(!BCP_GCG_low)

# Plot specification
plot_lm <- \(protein, centroid, label) {
	protein <- as.character(protein)

	y_max <- max(plot_data[[protein]], na.rm = T) * 1.25
	y_min <- min(plot_data[[protein]], na.rm = T)

	p <- ggplot(plot_data, aes(.data[[centroid]], .data[[protein]])) +
		facet_wrap("donor") +
		geom_point() +
		geom_smooth(method = "lm") +
		theme_bw(12) +
		ggpubr::stat_regline_equation() +
		scale_y_continuous(limits = c(y_min, y_max)) +
		labs(x = label, y = protein)

	png(
		str_glue("{base_path}/RD11a_1-{protein}_vs_{label}.png"),
		height = 6,
		width = 8,
		units = "in",
		res = 600
	)
	plot(p)
	dev.off()
}

# Hyaluronan-related proteins vs IIRS -------------------------------------

hyaluronan_specs <- tribble(
	~protein , ~centroid         , ~label ,
	"ITIH5"  , "im_sig_centroid" , "IIRS" ,
	"HEXA"   , "im_sig_centroid" , "IIRS" ,
	"HEXB"   , "im_sig_centroid" , "IIRS"
)

pwalk(hyaluronan_specs, plot_lm)


# Extracellular matrix (ECM) proteins vs BCP ------------------------------

ecm_specs <- tribble(
	~protein , ~centroid      , ~label ,
	"FAP"    , "bcp_centroid" , "BCP"  ,
	"DPP4"   , "bcp_centroid" , "BCP"  ,
	"QSOX1"  , "bcp_centroid" , "BCP"  ,
	"FBLN7"  , "bcp_centroid" , "BCP"
)

pwalk(ecm_specs, plot_lm)
