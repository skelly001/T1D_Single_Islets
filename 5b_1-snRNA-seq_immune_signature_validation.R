# ==============================================================================
# Immune Signature Validation Using snRNA-seq Reference Data
# ==============================================================================
# Script: 5b_1-snRNA-seq_immune_signature_validation.R
# Description: Validates the single-islet proteomics immune signature against
#              published snRNA-seq data from recent-onset T1D patients
#              (DOI:10.1126/sciadv.ady0080). Performs donor-level linear modeling
#              with limma and generates heatmap comparing log fold-changes between
#              proteomics (Stage 1 T1D) and reference transcriptomics (recent-onset T1D).
#
# Input: - output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_feature_candidates.xlsx
#        - data/Kyle_Gaulton_10.1126sciadv.ady0080_Beta_Cell_differentiations.xlsx
#        - output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds
# Output: - output/RD5-islet_immune_response_signature_analysis/RD5b_1-immune_signature_validation.png
# ==============================================================================

library(MSnSet.utils)
library(tidyverse)
library(ComplexHeatmap)


# Import ------------------------------------------------------------------

# Load immune signature candidate proteins
iirs_prot <- readxl::read_xlsx(
	"output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_feature_candidates.xlsx"
)

# Load published snRNA-seq data (DOI:10.1126/sciadv.ady0080)
ref_expr <- readxl::read_xlsx(
	"data/Kyle_Gaulton_10.1126sciadv.ady0080_Beta_Cell_differentiations.xlsx",
	sheet = "Early T1D_Significants"
)

# Load protein abundance with IIRS centroid
m <- readRDS("output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds")


# Process -----------------------------------------------------------------

## Donor-level IIRS correlation analysis
m_donors <- split(m, "donor")

iirs_lm_res <- vector("list", length(m_donors))
names(iirs_lm_res) <- names(m_donors)

for (i in names(m_donors)) {
	m_lp <- m_donors[[i]]

	# Test protein-IIRS associations via limma
	iirs_lm_res[[i]] <- limma_gen(m_lp, "~ im_sig_centroid", "im_sig_centroid") %>%
		rownames_to_column("gene_name") %>%
		arrange(P.Value) %>%
		left_join(fData(m_lp)[, c("gene_name", "ENTREZID")], by = "gene_name") %>%
		relocate(ENTREZID, .after = gene_name)
}


# Filter reference data for significant beta cell genes in early T1D
ref_expr_sig <- ref_expr %>%
	select(-`...1`) %>%
	filter(padj < 0.1, celltype == "Beta", state == "T1D_early")

# Identify genes present in both IIRS candidates and reference
gene_overlap <- ref_expr_sig[ref_expr_sig$Gene %in% iirs_prot$gene_name, ]

ref_logfc <- gene_overlap %>%
	select(Gene, log2FoldChange) %>%
	rename(Reference = log2FoldChange, gene_name = Gene)


# Prepare case donor logFC for comparison
prep_stats <- \(x, idx) {
	x %>%
		select(gene_name, logFC) %>%
		rename(!!str_glue("{idx}") := logFC)
}

case_donors <- c("6450", "6521", "6267")

iirs_lm_res_wide <- iirs_lm_res[case_donors] %>%
	imap(prep_stats) %>%
	reduce(left_join)

# Combine proteomics and reference transcriptomics data
plot_data <- inner_join(iirs_lm_res_wide, ref_logfc, "gene_name") %>%
	column_to_rownames("gene_name") %>%
	as.matrix()


## Visualization ----

# Color scale
col_pal <- circlize::colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

# Column labels
proteomics_label <- "Single-Islet\nProteomics\n(Stage 1 T1D)"
transcriptomics_label <- "snRNA-seq\n(Recent-onset T1D)"
col_order <- factor(
	c(rep(proteomics_label, 3), transcriptomics_label),
	levels = c(proteomics_label, transcriptomics_label)
)

# Remove column name for reference column
colnames(plot_data)[4] <- ""

hmap <- Heatmap(
	plot_data,
	column_split = col_order,
	column_gap = unit(20, "mm"),
	cluster_columns = F,
	column_names_rot = 0,
	column_names_centered = TRUE,
	show_row_dend = F,
	col = col_pal,
	name = "logFC",
	width = unit(16, "mm") * ncol(plot_data),
	height = unit(5, "mm") * nrow(plot_data)
)

png(
	"output/RD5-islet_immune_response_signature_analysis/RD5b_1-immune_signature_validation.png",
	width = 4.3,
	height = 3.6,
	units = "in",
	res = 600
)
draw(hmap)
dev.off()
