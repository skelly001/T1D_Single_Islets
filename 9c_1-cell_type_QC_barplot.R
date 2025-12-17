# ==============================================================================
# Cell Type Marker Quality Control Visualization
# ==============================================================================
# Script: 9c_1-cell_type_QC_barplot.R
# Description: Visualizes relative abundance of cell type markers across donors
#              using iBAQ (intensity-based absolute quantification) intensities.
#              Back-transforms log2 intensities, calculates relative expression by
#              donor, and generates horizontal stacked barplots showing cell type
#              composition. Creates composite figure with grid arrangement.
#
# Input: - output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_2-msnset_all_donors_BC_iBAQ.rds
#        - output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds
#        - output/RD9-cell_type_marker_QC/RD9b_1-scRNA-Seq_cell_type_markers.rds
# Output: - output/RD9-cell_type_marker_QC/RD9c_1-cell_type_marker_QC_barplot.png
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)
library(rcartocolor)
library(ggbreak)
library(gridExtra)
library(magick)


# Import ------------------------------------------------------------------

source("data/single_islets_project_functions.R")

# Load iBAQ protein abundance data
m <- readRDS(
	"output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_2-msnset_all_donors_BC_iBAQ.rds"
)

# Load IIRS protein features
im_sig <- readRDS("output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds")

# Load cell type markers from scRNA-seq
ctm <- readRDS("output/RD9-cell_type_marker_QC/RD9b_1-scRNA-Seq_cell_type_markers.rds")


# Process -----------------------------------------------------------------

# Calculate IIRS centroid
m_im <- m[fData(m)$UniProtAcc %in% im_sig$UniProtAcc, ]
m$im_sig_centroid <- apply(exprs(m_im), 2, \(x) {
	median(x, na.rm = T)
})

# Require 50% completeness per donor
m2 <- completeness_filt(m, 0.5)


## Prepare Cell Type Markers ----
ctm <- ctm %>%
	arrange(`Cell Type`, desc(Ratio)) %>%
	mutate(
		`Cell Type` = factor(
			`Cell Type`,
			levels = c(
				"alpha",
				"beta",
				"acinar",
				"delta",
				"gamma",
				"endothelial",
				"mast",
				"ductal",
				"epsilon",
				"activated_stellate",
				"macrophage"
			)
		)
	) %>%
	separate_wider_delim(rowname, delim = "_", names = c("gene_name", "UniProtAcc")) %>%
	select(-gene_name)

m3 <- m2[na.omit(match(ctm$UniProtAcc, fData(m2)$UniProtAcc)), ]

# Add cell type annotations to feature data
fData(m3) <- fData(m3) %>%
	rownames_to_column() %>%
	left_join(ctm, by = "UniProtAcc") %>%
	column_to_rownames()

# Select top 2 markers per cell type by mean intensity
sampnames <- paste0("X", sampleNames(m3))
top_marker_exp <- ms2df(m3) %>%
	rownames_to_column() %>%
	select(rowname, `Cell Type`, all_of(sampnames)) %>%
	rowwise() %>%
	mutate(mean_int = mean(c_across(where(is.numeric)), na.rm = T)) %>%
	ungroup() %>%
	slice_max(order_by = mean_int, n = 2, by = "Cell Type")

top_markers <- top_marker_exp$rowname

# Subset to top markers
m_ctm <- m3[fData(m3)$gene_name %in% top_markers, ]

# Back-transform from log2 scale
exprs(m_ctm) <- 2^exprs(m_ctm)

m_split <- split(m_ctm, "donor")

marker_mean_fun <- \(x) {
	ms2df(x) %>%
		select(gene_name, `Cell Type`, all_of(paste0("X", sampleNames(x)))) %>%
		rowwise() %>%
		mutate(mean_exp = mean(c_across(where(is.numeric)), na.rm = T)) %>%
		select(gene_name, `Cell Type`, mean_exp) %>%
		ungroup()
}

# Calculate mean expression per donor
marker_means <- map(m_split, marker_mean_fun) %>%
	enframe("donor") %>%
	unnest(value)

# Calculate relative expression per donor
marker_relative <- marker_means %>%
	mutate(donor_sum = sum(mean_exp, na.rm = T), .by = "donor") %>%
	mutate(relative_exp = mean_exp / donor_sum) %>%
	group_by(donor) %>%
	arrange(`Cell Type`, .by_group = T) %>%
	ungroup() %>%
	mutate(
		gene_name = fct(gene_name),
		donor = fct(donor, levels = c("6450", "6521", "6267", "6178", "6440", "6539"))
	)


## Visualization ----------------------------------------------------------

# Define y-axis label order
label_order <- marker_relative %>%
	filter(donor == "6450") %>%
	arrange(fct_rev(gene_name)) %>%
	pull(`Cell Type`)

# Define cell type colors
safe_pal_full <- carto_pal(11, "Safe")
names(safe_pal_full) <- levels(marker_relative$`Cell Type`)
safe_pal2 <- safe_pal_full[label_order]


donors <- levels(marker_relative$donor)

# Helper function: create barplot for a single donor with conditional formatting
plot_fun <- function(donor) {
	df <- marker_relative[marker_relative$donor == donor, ]

	legend_pos <- if (donor == "6539") "right" else "none"
	y_text <- if (donor == "6450") element_text(color = safe_pal2) else element_blank()
	y_ticks <- if (donor == "6450") element_line() else element_blank()
	y_label <- if (donor == "6450") "Cell Type Marker" else NULL
	x_label <- NULL

	ggplot(df, aes(x = relative_exp, y = fct_rev(gene_name), fill = `Cell Type`)) +
		geom_bar(stat = "identity") +
		scale_x_break(c(0.15, 0.65), expand = FALSE, space = 0.1) +
		scale_fill_manual(values = safe_pal_full) +
		scale_x_continuous(breaks = seq(0, 1, 0.1)) +
		theme_classic() +
		theme(
			axis.text.x.top = element_blank(),
			axis.ticks.x.top = element_blank(),
			axis.line.x.top = element_blank(),
			axis.text.y = y_text,
			axis.ticks.y = y_ticks,
			panel.spacing = unit(1, "lines"),
			plot.margin = unit(c(0, 0, 0, 0), "pt"),
			legend.position = legend_pos,
			panel.grid.major = element_line(color = "grey90")
		) +
		ylab(y_label) +
		xlab(x_label)
}

# Generate individual donor plots
plots <- map(donors, plot_fun)

widths <- c(2.2, 1.5, 1.5, 1.5, 1.5, 3.1)
iwalk(plots, \(x, idx) {
	width <- widths[idx]
	ggsave(
		str_glue("output/RD9-cell_type_marker_QC/RD9c_1-cell_type_QC_plot_{idx}_temp.png"),
		plot = x,
		width = width,
		height = 4,
		dpi = 600
	)
})


## Combine Plots ----
img_files <- list.files("output/RD9-cell_type_marker_QC", "RD9c_1", full.names = T)
imgs <- map(img_files, image_read)

# Combine horizontally
combined_img <- image_append(do.call(c, imgs), stack = FALSE)

image_write(combined_img, "output/RD9-cell_type_marker_QC/RD9c_1-cell_type_marker_QC_barplot.png")

# Clean up temporary files
file.remove(img_files)
