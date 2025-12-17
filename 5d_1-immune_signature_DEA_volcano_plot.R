# ==============================================================================
# Immune Signature Differential Expression Volcano Plots
# ==============================================================================
# Script: 5d_1-immune_signature_DEA_volcano_plot.R
# Description: Visualizes differential expression results from immune signature
#              analysis as volcano plots with significance thresholds (BH-adjusted
#              p-value < 0.01). Creates multi-panel faceted plots showing log2
#              fold-change versus -log10(adjusted p-value) for each donor with
#              proteins classified as positive, negative, or not significant.
#
# Input: - output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv
# Output: - output/RD5-islet_immune_response_signature_analysis/RD5d_1-immune_signature_DEA_volcano_plots.png
# ==============================================================================

library(TMSig)
library(MSnSet.utils)
library(tidyverse)


# Import -----------------------------------------------------------------

# Load limma differential expression results
limma_res <- read.csv(
	"output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv"
)


# Process ----------------------------------------------------------------

## Prepare Data ----
x <- limma_res %>%
	pivot_longer(
		cols = -protein,
		names_to = c(".value", "donor"),
		names_pattern = "(.*)_(.*)"
	) %>%
	mutate(
		donor = factor(donor, levels = unique(donor)),
		z = sign(logFC) * qnorm(P.Value / 2, lower.tail = FALSE),
		pt_color = case_when(
			adj.P.Val >= 0.01 ~ "NS",
			logFC < 0 ~ "Negative",
			TRUE ~ "Positive"
		),
		pt_color = factor(pt_color, levels = c("NS", "Negative", "Positive"))
	) %>%
	arrange(donor, pt_color)

# Summarize significant protein counts per donor
label_df <- x %>%
	summarise(
		.by = donor,
		n_pos = sum(pt_color == "Positive"),
		n_neg = sum(pt_color == "Negative"),
		n_no_change = sum(pt_color == "NS")
	) %>%
	mutate(
		label = sprintf(
			"Negative: %d\nPositive: %d",
			n_neg,
			n_pos
		),
		x = -1,
		y = 45
	)


## Volcano Plot ----
p <- ggplot(x, aes(x = logFC, y = -log10(adj.P.Val))) +
	geom_point(aes(color = pt_color), size = 1.2, shape = 16) +
	geom_label(
		aes(x = x, y = y, label = label),
		data = label_df,
		lineheight = 1,
		size = 8,
		size.unit = "pt",
		label.padding = unit(0.16, "lines"),
		hjust = 0,
		vjust = 1
	) +
	facet_wrap(~donor, ncol = 3L) +
	scale_color_manual(
		name = "Slope:",
		values = c("blue", "red", "grey90"),
		breaks = c("Negative", "Positive", "NS"),
		limits = c("Negative", "Positive")
	) +
	labs(
		x = latex2exp::TeX(
			"log$_2$(fold-change)"
		),
		y = latex2exp::TeX(
			"-log$_{10}$(BH Adjusted P-value)"
		)
	) +
	theme_bw() +
	theme(
		panel.grid.minor.y = element_blank(),
		legend.position = "top",
		legend.direction = "horizontal",
		legend.margin = margin(b = -6)
	)


png(
	filename = "output/RD5-islet_immune_response_signature_analysis/RD5d_1-immune_signature_DEA_volcano_plots.png",
	height = 4,
	width = 5.5,
	units = "in",
	res = 600
)
p
dev.off()
