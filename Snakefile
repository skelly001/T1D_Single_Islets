# Snakemake workflow for T1D Single Islets proteomics analysis


# Final target - key outputs from all analysis stages
rule all:
	input:
		# Ensure renv environment is set up first
		"renv/.restored",
		# Stage 1 - MSnSets
		"output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData",
		"output/RD1-raw_msnsets/RD1a_1-msnsets_ref_raw.RData",
		# Stage 2 - preprocessing
		"output/RD2-preprocessing/RD2a_1-case_donors_preprocessing/RD2a_1-msnsets_donor_1_2_3_preprocessed.rds",
		"output/RD2-preprocessing/RD2a_2-control_donors_preprocessing/RD2a_2-msnsets_donor_4_5_6_preprocessed.rds",
		"output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds",
		"output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds",
		"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-msnsets_donor_1_2_3_preprocessed.rds",
		"output/RD2-preprocessing/RD2c_2-control_donors_preprocessing_iBAQ/RD2c_2-msnsets_donor_4_5_6_preprocessed.rds",
		"output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_1-msnset_all_donors_iBAQ.rds",
		"output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_2-msnset_all_donors_BC_iBAQ.rds",
		# Stage 3 - WGCNA
		"output/RD3-WGCNA/donor_6450/SoftPower_plot.png",
		"output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_module_membership.xlsx",
		"output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_ORA_BP_results.xlsx",
		"output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_ORA_dotplot.png",
		"output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_module_correlation_heatmap.png",
		"output/RD3-WGCNA/donor_6521/SoftPower_plot.png",
		"output/RD3-WGCNA/donor_6521/RD3b_1-WGCNA_module_membership.xlsx",
		"output/RD3-WGCNA/donor_6521/RD3b_1-WGCNA_ORA_BP_results.xlsx",
		"output/RD3-WGCNA/donor_6521/RD3b_1-WGCNA_ORA_dotplot.png",
		"output/RD3-WGCNA/donor_6521/RD3b_1-WGCNA_module_correlation_heatmap.png",
		"output/RD3-WGCNA/donor_6267/SoftPower_plot.png",
		"output/RD3-WGCNA/donor_6267/RD3c_1-WGCNA_module_membership.xlsx",
		"output/RD3-WGCNA/donor_6267/RD3c_1-WGCNA_ORA_BP_results.xlsx",
		"output/RD3-WGCNA/donor_6267/RD3c_1-WGCNA_ORA_dotplot.png",
		"output/RD3-WGCNA/donor_6267/RD3c_1-WGCNA_module_correlation_heatmap.png",
		# Stage 4 - immune signature
		"output/RD4-islet_immune_response_signature/RD4a_1-WGCNA_module_immune_signature_candidate_proteins.xlsx",
		"output/RD4-islet_immune_response_signature/RD4b_1-final_ML_models.rds",
		"output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds",
		"output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_feature_candidates.xlsx",
		"output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds",
		"output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds",
		"output/RD4-islet_immune_response_signature/RD4c_1-islet_immune_response_signature_heatmap.png",
		# Stage 5 - immune signature analysis
		"output/RD5-islet_immune_response_signature_analysis/RD5a_1-UMAP_IIRS_trajectory.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5a_1-UMAP_IIRS_trajectory_donor_faceted.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5b_1-immune_signature_validation.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv",
		"output/RD5-islet_immune_response_signature_analysis/RD5d_1-immune_signature_DEA_volcano_plots.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5e_1-gene_set_table.csv",
		"output/RD5-islet_immune_response_signature_analysis/RD5e_1-immune_signature_CAMERA-PR_significant_results.csv",
		"output/RD5-islet_immune_response_signature_analysis/RD5e_1-immune_signature_CAMERA-PR_bubble_plot.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5f_1-ECM-related_GO_term_heatmap.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5f_1-hyaluronan_metabolic_process_GO_term_heatmap.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5f_1-glycosaminoglycan_catabolic_process_GO_term_heatmap.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5f_1-interleukin-10_production_GO_term_heatmap.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_intradonor_variance_boxplot.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_STAT1_WARS1_correlation_scatterplots.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_correlation_heatmap.png",
		# Stage 6 - beta cell profile
		"output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx",
		"output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds",
		"output/RD6-beta_cell_profile/RD6b_1-beta_cell_profile_heatmap.png",
		# Stage 7 - beta cell profile analysis
		"output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv",
		"output/RD7-beta_cell_profile_analysis/RD7b_1-beta_cell_profile_DEA_volcano_plots.png",
		"output/RD7-beta_cell_profile_analysis/RD7c_1-gene_set_table.csv",
		"output/RD7-beta_cell_profile_analysis/RD7c_1-beta_cell_profile_CAMERA-PR_significant_results.csv",
		"output/RD7-beta_cell_profile_analysis/RD7c_1-beta_cell_profile_CAMERA-PR_bubble_plot.png",
		"output/RD7-beta_cell_profile_analysis/RD7d_1-UMAP_IIRS_BCP_trajectory_comparison_donor_faceted.png",
		"output/RD7-beta_cell_profile_analysis/RD7e_1-BCP_intradonor_variance_boxplot.png",
		# Stage 8 - clustering
		"output/RD8-clustering/RD8b_1-clustering_results.RData",
		"output/RD8-clustering/RD8b_1-cluster_membership.xlsx",
		"output/RD8-clustering/RD8c_1-cluster_heatmaps.pdf",
		"output/RD8-clustering/RD8c_1-LGALS3_LGALS3BP_heatmap.pdf",
		"output/RD8-clustering/RD8d_1-basement_membrane_vs_BCP.png",
		"output/RD8-clustering/RD8d_1-basement_membrane_vs_IIRS.png",
		# Stage 10 - final QC
		"output/RD10-misc/RD10a_1-INS_intensity_boxplot.png",
		"output/RD10-misc/RD10a_1-GCG_intensity_boxplot.png",
		"output/RD10-misc/RD10a_1-INS_vs_GCG_scatterplot.png",
		"output/RD10-misc/RD10a_1-observed_proteins_QC_barplot.png",
		"output/RD10-misc/RD10b_1-table_1_donor_information.png",
		# Stage 11 - protein-centroid linear modeling
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-ITIH5_vs_IIRS.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-HEXA_vs_IIRS.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-HEXB_vs_IIRS.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-FAP_vs_BCP.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-DPP4_vs_BCP.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-QSOX1_vs_BCP.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-FBLN7_vs_BCP.png"

# ==============================================================================
# Setup: Initialize renv environment
# ==============================================================================

rule setup_renv:
	input:
		"renv.lock",
		".Rprofile"
	output:
		touch("renv/.restored")
	shell:
		'Rscript -e "renv::restore(prompt = FALSE)"'

# ==============================================================================
# Stage 1: Create MSnSets from raw MSstats data
# ==============================================================================

rule create_msnsets:
	input:
		script = "1a_1-create_msnsets.R",
		metadata = "data/sample_metadata.xlsx"
	output:
		"output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData",
		"output/RD1-raw_msnsets/RD1a_1-msnsets_ref_raw.RData"
	shell:
		"Rscript {input.script}"

# ==============================================================================
# Stage 2: Preprocessing - donor-level processing and combination
# ==============================================================================

# Process case donors (1-3: T1D donors 6450, 6521, 6267)
rule process_msnsets_donors_1_2_3:
	input:
		script = "2a_1-process_msnsets_donors_1_2_3.R",
		function_script = "data/single_islets_project_functions.R",
		msnsets = "output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData",
		uniprot = "data/uniprot_database_observed_proteins.tsv"
	output:
		"output/RD2-preprocessing/RD2a_1-case_donors_preprocessing/RD2a_1-msnsets_donor_1_2_3_preprocessed.rds"
	shell:
		"Rscript {input.script}"

# Process control donors (4-6: non-diabetic donors 6178, 6440, 6539)
rule process_msnsets_donors_4_5_6:
	input:
		script = "2a_2-process_msnsets_donors_4_5_6.R",
		function_script = "data/single_islets_project_functions.R",
		msnsets = "output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData",
		uniprot = "data/uniprot_database_observed_proteins.tsv"
	output:
		"output/RD2-preprocessing/RD2a_2-control_donors_preprocessing/RD2a_2-msnsets_donor_4_5_6_preprocessed.rds"
	shell:
		"Rscript {input.script}"

# Combine all donors without donor-level batch correction
rule combine_msnsets_all_donors:
	input:
		script = "2b_1-combine_msnsets_all_donors.R",
		case = "output/RD2-preprocessing/RD2a_1-case_donors_preprocessing/RD2a_1-msnsets_donor_1_2_3_preprocessed.rds",
		control = "output/RD2-preprocessing/RD2a_2-control_donors_preprocessing/RD2a_2-msnsets_donor_4_5_6_preprocessed.rds"
	output:
		"output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds"
	shell:
		"Rscript {input.script}"

# Combine all donors with donor-level batch correction
rule combine_msnsets_batch_correction:
	input:
		script = "2b_2-combine_msnsets_all_donors_with_donor_batch_correction.R",
		case = "output/RD2-preprocessing/RD2a_1-case_donors_preprocessing/RD2a_1-msnsets_donor_1_2_3_preprocessed.rds",
		control = "output/RD2-preprocessing/RD2a_2-control_donors_preprocessing/RD2a_2-msnsets_donor_4_5_6_preprocessed.rds"
	output:
		"output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds"
	shell:
		"Rscript {input.script}"

# Process case donors with iBAQ normalization
rule process_msnsets_donors_1_2_3_iBAQ:
	input:
		script = "2c_1-process_msnsets_donors_1_2_3_iBAQ.R",
		function_script = "data/single_islets_project_functions.R",
		msnsets = "output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData",
		uniprot = "data/uniprot_database_observed_proteins.tsv",
		fasta = "data/Results_7_combined/2024-10-09-decoys-contam-uniprotkb_Human_AND_reviewed_true_AND_m_2024_10_09.fasta.fas"
	output:
		"output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-msnsets_donor_1_2_3_preprocessed.rds"
	shell:
		"Rscript {input.script}"

# Process control donors with iBAQ normalization
rule process_msnsets_donors_4_5_6_iBAQ:
	input:
		script = "2c_2-process_msnsets_donors_4_5_6_iBAQ.R",
		function_script = "data/single_islets_project_functions.R",
		msnsets = "output/RD1-raw_msnsets/RD1a_1-msnsets_raw.RData",
		uniprot = "data/uniprot_database_observed_proteins.tsv",
		fasta = "data/Results_7_combined/2024-10-09-decoys-contam-uniprotkb_Human_AND_reviewed_true_AND_m_2024_10_09.fasta.fas"
	output:
		"output/RD2-preprocessing/RD2c_2-control_donors_preprocessing_iBAQ/RD2c_2-msnsets_donor_4_5_6_preprocessed.rds"
	shell:
		"Rscript {input.script}"

# Combine all donors with iBAQ, no donor batch correction
rule combine_msnsets_all_donors_iBAQ:
	input:
		script = "2d_1-combine_msnsets_all_donors_iBAQ.R",
		case = "output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-msnsets_donor_1_2_3_preprocessed.rds",
		control = "output/RD2-preprocessing/RD2c_2-control_donors_preprocessing_iBAQ/RD2c_2-msnsets_donor_4_5_6_preprocessed.rds"
	output:
		"output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_1-msnset_all_donors_iBAQ.rds"
	shell:
		"Rscript {input.script}"

# Combine all donors with iBAQ and donor batch correction
rule combine_msnsets_batch_correction_iBAQ:
	input: 
		script = "2d_2-combine_msnsets_all_donors_with_donor_batch_correction_iBAQ.R",
		case = "output/RD2-preprocessing/RD2c_1-case_donors_preprocessing_iBAQ/RD2c_1-msnsets_donor_1_2_3_preprocessed.rds",
		control = "output/RD2-preprocessing/RD2c_2-control_donors_preprocessing_iBAQ/RD2c_2-msnsets_donor_4_5_6_preprocessed.rds"
	output:
		"output/RD2-preprocessing/RD2d-prepared_data_iBAQ/RD2d_2-msnset_all_donors_BC_iBAQ.rds"
	shell:
		"Rscript {input.script}"

# ==============================================================================
# Stage 3: WGCNA co-expression network analysis (one per mAAb+ donor)
# ==============================================================================

rule wgcna_donor_6450:
	input:
		script = "3a_1-WGCNA_donor_6450.R",
		function_script = "data/GSEA_functions.R",
		msnset = "output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds"
	output:
		"output/RD3-WGCNA/donor_6450/SoftPower_plot.png",
		"output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_module_membership.xlsx",
		"output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_ORA_BP_results.xlsx",
		"output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_ORA_dotplot.png",
		"output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_module_correlation_heatmap.png"
	shell:
		"Rscript {input.script}"

rule wgcna_donor_6521:
	input:
		script = "3b_1-WGCNA_donor_6521.R",
		function_script = "data/GSEA_functions.R",
		msnset = "output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds"
	output:
		"output/RD3-WGCNA/donor_6521/SoftPower_plot.png",
		"output/RD3-WGCNA/donor_6521/RD3b_1-WGCNA_module_membership.xlsx",
		"output/RD3-WGCNA/donor_6521/RD3b_1-WGCNA_ORA_BP_results.xlsx",
		"output/RD3-WGCNA/donor_6521/RD3b_1-WGCNA_ORA_dotplot.png",
		"output/RD3-WGCNA/donor_6521/RD3b_1-WGCNA_module_correlation_heatmap.png"
	shell:
		"Rscript {input.script}"

rule wgcna_donor_6267:
	input:
		script = "3c_1-WGCNA_donor_6267.R",
		function_script = "data/GSEA_functions.R",
		msnset = "output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds"
	output:
		"output/RD3-WGCNA/donor_6267/SoftPower_plot.png",
		"output/RD3-WGCNA/donor_6267/RD3c_1-WGCNA_module_membership.xlsx",
		"output/RD3-WGCNA/donor_6267/RD3c_1-WGCNA_ORA_BP_results.xlsx",
		"output/RD3-WGCNA/donor_6267/RD3c_1-WGCNA_ORA_dotplot.png",
		"output/RD3-WGCNA/donor_6267/RD3c_1-WGCNA_module_correlation_heatmap.png"
	shell:
		"Rscript {input.script}"

# ==============================================================================
# Stage 4: Islet Immune Response Signature (IIRS) identification
# ==============================================================================

rule immune_signature_candidate_proteins:
	input:
		script = "4a_1-immune_signature_candidate_proteins.R",
		wgcna_6450 = "output/RD3-WGCNA/donor_6450/RD3a_1-WGCNA_module_membership.xlsx",
		wgcna_6521 = "output/RD3-WGCNA/donor_6521/RD3b_1-WGCNA_module_membership.xlsx",
		wgcna_6267 = "output/RD3-WGCNA/donor_6267/RD3c_1-WGCNA_module_membership.xlsx",
		msnset = "output/RD2-preprocessing/RD2b-prepared_data/RD2b_1-msnset_all_donors.rds"
	output:
		"output/RD4-islet_immune_response_signature/RD4a_1-WGCNA_module_immune_signature_candidate_proteins.xlsx"
	shell:
		"Rscript {input.script}"

rule immune_signature_selection:
	input:
		script = "4b_1-immune_signature_selection.R",
		function_script = "data/single_islets_project_functions.R",
		msnset = "output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds",
		candidates = "output/RD4-islet_immune_response_signature/RD4a_1-WGCNA_module_immune_signature_candidate_proteins.xlsx"
	output:
		"output/RD4-islet_immune_response_signature/RD4b_1-final_ML_models.rds",
		"output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds",
		"output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_feature_candidates.xlsx",
		"output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds"
	shell:
		"Rscript {input.script}"

rule immune_signature_heatmap:
	input:
		script = "4c_1-immune_signature_heatmap.R",
		function_script = "data/single_islets_project_functions.R",
		msnset = "output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds",
		features = "output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds"
	output:
		"output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds",
		"output/RD4-islet_immune_response_signature/RD4c_1-islet_immune_response_signature_heatmap.png"
	shell:
		"Rscript {input.script}"

# ==============================================================================
# Stage 5: Immune signature analysis and validation
# ==============================================================================

rule immune_signature_trajectory_analysis:
	input:
		script = "5a_1-immune_signature_trajectory_analysis.R",
		function_script = "data/single_islets_project_functions.R",
		msnset_bc = "output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds",
		msnset_im = "output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds"
	output:
		"output/RD5-islet_immune_response_signature_analysis/RD5a_1-UMAP_IIRS_trajectory.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5a_1-UMAP_IIRS_trajectory_donor_faceted.png"
	shell:
		"Rscript {input.script}"

rule snRNA_seq_immune_signature_validation:
	input:
		script = "5b_1-snRNA-seq_immune_signature_validation.R",
		candidates = "output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_feature_candidates.xlsx",
		msnset = "output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds",
		RNAseq = "data/Kyle_Gaulton_10.1126sciadv.ady0080_Beta_Cell_differentiations.xlsx"
	output:
		"output/RD5-islet_immune_response_signature_analysis/RD5b_1-immune_signature_validation.png"
	shell:
		"Rscript {input.script}"

rule immune_signature_differential_analysis:
	input:
		script = "5c_1-immune_signature_differential_analysis.R",
		function_script = "data/single_islets_project_functions.R",
		msnset = "output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds"
	output:
		"output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv"
	shell:
		"Rscript {input.script}"

rule immune_signature_DEA_volcano_plot:
	input:
		script = "5d_1-immune_signature_DEA_volcano_plot.R",
		limma = "output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv"
	output:
		"output/RD5-islet_immune_response_signature_analysis/RD5d_1-immune_signature_DEA_volcano_plots.png"
	shell:
		"Rscript {input.script}"

rule immune_signature_CAMERA_PR:
	input:
		script = "5e_1-immune_signature_CAMERA-PR.R",
		limma = "output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv",
		gene_sets = "data/GOBP_gene_sets.rds"
	output:
		"output/RD5-islet_immune_response_signature_analysis/RD5e_1-gene_set_table.csv",
		"output/RD5-islet_immune_response_signature_analysis/RD5e_1-immune_signature_CAMERA-PR_significant_results.csv",
		"output/RD5-islet_immune_response_signature_analysis/RD5e_1-immune_signature_CAMERA-PR_bubble_plot.png"
	shell:
		"Rscript {input.script}"

rule GO_term_heatmaps:
	input:
		script = "5f_1-GO_term_heatmaps.R",
		limma = "output/RD5-islet_immune_response_signature_analysis/RD5c_1-immune_signature_limma_results.csv",
		gene_sets = "data/GOBP_gene_sets.rds"
	output:
		"output/RD5-islet_immune_response_signature_analysis/RD5f_1-ECM-related_GO_term_heatmap.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5f_1-hyaluronan_metabolic_process_GO_term_heatmap.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5f_1-glycosaminoglycan_catabolic_process_GO_term_heatmap.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5f_1-interleukin-10_production_GO_term_heatmap.png"
	shell:
		"Rscript {input.script}"

rule immune_signature_misc_plots:
	input:
		script = "5g_1-immune_signature_misc_plots.R",
		msnset = "output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds",
		msnset_im = "output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds",
		features = "output/RD4-islet_immune_response_signature/RD4b_1-immune_signature_features.rds"
	output:
		"output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_intradonor_variance_boxplot.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_STAT1_WARS1_correlation_scatterplots.png",
		"output/RD5-islet_immune_response_signature_analysis/RD5g_1-IIRS_correlation_heatmap.png"
	shell:
		"Rscript {input.script}"

# ==============================================================================
# Stage 6: Beta cell profile identification
# ==============================================================================

rule beta_cell_profile_selection:
	input:
		script = "6a_1-beta_cell_profile_selection.R",
		function_script = "data/single_islets_project_functions.R",
		msnset = "output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds",
		msnset_im = "output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds"
	output:
		"output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx",
		"output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds"
	shell:
		"Rscript {input.script}"

rule beta_cell_profile_heatmap:
	input:
		script = "6b_1-beta_cell_profile_heatmap.R",
		function_script = "data/single_islets_project_functions.R",
		msnset = "output/RD4-islet_immune_response_signature/RD4b_1-msnset_all_donors_im_sig.rds",
		msnset_im = "output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds",
		features = "output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx"
	output:
		"output/RD6-beta_cell_profile/RD6b_1-beta_cell_profile_heatmap.png"
	shell:
		"Rscript {input.script}"

# ==============================================================================
# Stage 7: Beta cell profile analysis
# ==============================================================================

rule beta_cell_profile_differential_analysis:
	input:
		script = "7a_1-beta_cell_profile_differential_analysis.R",
		function_script = "data/single_islets_project_functions.R",
		msnset = "output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds"
	output:
		"output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv"
	shell:
		"Rscript {input.script}"

rule beta_cell_profile_DEA_volcano_plot:
	input:
		script = "7b_1-beta_cell_profile_DEA_volcano_plot.R",
		limma = "output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv"
	output:
		"output/RD7-beta_cell_profile_analysis/RD7b_1-beta_cell_profile_DEA_volcano_plots.png"
	shell:
		"Rscript {input.script}"

rule beta_cell_profile_CAMERA_PR:
	input:
		script = "7c_1-beta_cell_profile_CAMERA-PR.R",
		limma = "output/RD7-beta_cell_profile_analysis/RD7a_1-beta_cell_profile_limma_results.csv",
		gene_sets = "data/GOBP_gene_sets.rds"
	output:
		"output/RD7-beta_cell_profile_analysis/RD7c_1-gene_set_table.csv",
		"output/RD7-beta_cell_profile_analysis/RD7c_1-beta_cell_profile_CAMERA-PR_significant_results.csv",
		"output/RD7-beta_cell_profile_analysis/RD7c_1-beta_cell_profile_CAMERA-PR_bubble_plot.png"
	shell:
		"Rscript {input.script}"

rule beta_cell_profile_IIRS_trajectory_comparison:
	input:
		script = "7d_1-beta_cell_profile_IIRS_trajectory_comparison.R",
		function_script = "data/single_islets_project_functions.R",
		msnset = "output/RD2-preprocessing/RD2b-prepared_data/RD2b_2-msnset_all_donors_BC.rds",
		msnset_im = "output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds",
		features = "output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx"
	output:
		"output/RD7-beta_cell_profile_analysis/RD7d_1-UMAP_IIRS_BCP_trajectory_comparison_donor_faceted.png"
	shell:
		"Rscript {input.script}"

rule beta_cell_profile_misc_plots:
	input:
		script = "7e_1-beta_cell_profile_misc_plots.R",
		msnset = "output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds",
		features = "output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx"
	output:
		"output/RD7-beta_cell_profile_analysis/RD7e_1-BCP_intradonor_variance_boxplot.png"
	shell:
		"Rscript {input.script}"

# ==============================================================================
# Stage 8: Clustering analysis
# ==============================================================================

rule clustering:
	input:
		script = "8b_1-clustering.R",
		function_script = "8a_1-clustGO_Function_GOScore.R",
		msnset = "output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds"
	output:
		"output/RD8-clustering/RD8b_1-clustering_results.RData",
		"output/RD8-clustering/RD8b_1-cluster_membership.xlsx"
	shell:
		"Rscript {input.script}"

rule cluster_heatmaps:
	input:
		script = "8c_1-cluster_heatmaps.R",
		function_script = "data/single_islets_project_functions.R",
		clustering = "output/RD8-clustering/RD8b_1-clustering_results.RData",
		msnset = "output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds",
		msnset_im = "output/RD4-islet_immune_response_signature/RD4c_1-msnset_immune_sig_case_adjusted.rds",
		features = "output/RD6-beta_cell_profile/RD6a_1-beta_cell_profile_features.xlsx"
	output:
		"output/RD8-clustering/RD8c_1-cluster_heatmaps.pdf",
		"output/RD8-clustering/RD8c_1-LGALS3_LGALS3BP_heatmap.pdf"
	shell:
		"Rscript {input.script}"

rule basement_membrane_cluster:
	input:
		script = "8d_1-basement_membrane_cluster.R",
		clustering = "output/RD8-clustering/RD8b_1-clustering_results.RData",
		msnset = "output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds"
	output:
		"output/RD8-clustering/RD8d_1-basement_membrane_vs_BCP.png",
		"output/RD8-clustering/RD8d_1-basement_membrane_vs_IIRS.png"
	shell:
		"Rscript {input.script}"

# ==============================================================================
# Stage 9: Cell type marker quality control
# ==============================================================================

# 1. Run Azimuth human pancreas scRNA-seq snakefile separately.
# File Path: "./azimuth-references/human_pancreas_snakemake/Snakefile"
# 2. Run scripts 9a_1, 9b_1, 9c_1.

# ==============================================================================
# Stage 10: Final QC plots and summary tables
# ==============================================================================

rule misc_QC_plots:
	input:
		script = "10a_1-misc_QC_plots.R",
		msnset = "output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds"
	output:
		"output/RD10-misc/RD10a_1-INS_intensity_boxplot.png",
		"output/RD10-misc/RD10a_1-GCG_intensity_boxplot.png",
		"output/RD10-misc/RD10a_1-INS_vs_GCG_scatterplot.png",
		"output/RD10-misc/RD10a_1-observed_proteins_QC_barplot.png"
	shell:
		"Rscript {input.script}"

rule table_1_donor_information:
	input:
		script = "10b_1_table_1_donor_information.R",
		data = "data/donor_table.xlsx"
	output:
		"output/RD10-misc/RD10b_1-table_1_donor_information.png"
	shell:
		"Rscript {input.script}"

# ==============================================================================
# Stage 11: Protein–centroid linear modeling
# ==============================================================================

rule protein_centroid_linear_modeling:
	input:
		script = "11a_1-protein_centroid_linear_modeling.R",
		msnset = "output/RD6-beta_cell_profile/RD6a_1-msnset_w_centroids.rds"
	output:
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-ITIH5_vs_IIRS.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-HEXA_vs_IIRS.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-HEXB_vs_IIRS.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-FAP_vs_BCP.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-DPP4_vs_BCP.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-QSOX1_vs_BCP.png",
		"output/RD11-protein_centroid_linear_modeling/RD11a_1-FBLN7_vs_BCP.png"
	shell:
		"Rscript {input.script}"
