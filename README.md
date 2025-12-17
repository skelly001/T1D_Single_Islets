# Single-Islet Proteomics Maps Pseudo-Temporal Islet Immune Responses and Dysfunction in Stage 1 Type 1 Diabetes

[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.1101%2F2025.11.10.687674-b31b1b)](https://doi.org/10.1101/2025.11.10.687674) [![MassIVE](https://img.shields.io/badge/MassIVE-MSV000099920-blue)](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=MSV000099920)

**Shane S. Kelly, Soumyadeep Sarkar, Sarai M. Williams, An D. Fu, Elizabeth A. Butterworth, Tyler J. Sagendorf, Lorenz A. Nierves, Yumi Kwon, Xiaolu Li, Vladislav A. Petyuk, James M. Fulcher, Jing Chen, Ernesto S. Nakayasu, Mark A. Atkinson, Rohit N. Kulkarni, Clayton E Mathews, Ying Zhu, Martha Campbell-Thompson, and Wei-Jun Qian**

This repository contains the complete analysis pipeline for single-islet proteomics in stage 1 Type 1 Diabetes (T1D), supporting our manuscript currently available as a [preprint](https://www.biorxiv.org/content/10.1101/2025.11.10.687674v1.full).

## Overview

Progressive β-cell dysfunction precedes the onset of T1D, yet the molecular mechanisms driving early T1D development remain poorly understood. This study applies a single-islet proteomics workflow to profile intra-donor islet heterogeneity in three stage 1 T1D cases with matched non-diabetic controls, defining *in situ* protein signatures of pseudo-temporal islet dysfunction.

### Key Findings

-   **\~100 individual islets per donor** were analyzed using laser microdissection (LMD) and nanoPOTS-based proteomics, revealing consistent proteomic patterns reflecting pseudo-time progression.
-   **Islet Immune Response Signature (IIRS)**: A 40-protein signature capturing immune-mediated islet progression, including HLA class I/II machinery, interferon-stimulated genes, antigen processing components, and novel candidates (PARP10, GSDMD, LGALS3BP, OPTN).
-   **Beta Cell Profile (BCP)**: A 42-protein panel representing β-cell identity and function, revealing progressive loss of β-cell characteristics in stage 1 T1D.
-   **Extracellular matrix (ECM) dysregulation** was identified in association with islet immune response, while **mRNA processing and RNA-splicing** pathways were prominently associated with loss of β-cell function. ECM remodeling was also negatively correlated with β-cell function.
-   Evidence of pseudo-temporal trajectories within individual donors, with immune activation and β-cell dysfunction occurring largely independently (Pearson correlation r = 0.13).

### Study Design

-   **Stage 1 T1D Donors**: 3 donors (6450, 6521, 6267) - multiple autoantibody positive (mAAb+)
-   **Non-diabetic Controls**: 3 age, sex, and race-matched donors (6178, 6440, 6539)
-   **Methodology**: Single-islet spatial proteomics (\~5,800 proteins per donor), multiplex immunohistochemistry (mIHC), machine learning, WGCNA, and pathway enrichment analyses

## Data Availability

### Proteomic Data

Raw mass spectrometry data and processed MSstats files are available from MassIVE: - **Accession**: [MSV000099920](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=MSV000099920)

Required files from MassIVE: - `MSstats.csv` output files (FragPipe processed quantification data) - `sample_metadata.xlsx` (sample phenotype and experimental metadata) - `uniprotkb_Human_2023_10_25.tsv.gz` (UniProt human proteome reference for protein annotation)

### Data Setup

Before running the analysis pipeline, download the required files from MassIVE and place them in the appropriate locations:

``` bash
# Create data directory structure
mkdir data/Results_7_combined

# Download all files from MassIVE and place them in the data directories
# - MSstats.csv files → data/Results_7_combined/Results_7.1/, Results_7.2/, and Results_7.3/
# - sample_metadata.xlsx → data/
# - uniprotkb_Human_2023_10_25.tsv.gz → data/
```

Additional required reference files are included in this repository under the data folder.

## Installation

### Prerequisites

-   R (≥4.5)
-   Snakemake (≥9.14)

### Environment Setup

This project uses `renv` for R package management, ensuring reproducible package versions. The Snakefile automatically initializes the environment on first run.

``` bash
# Clone the repository
git clone https://github.com/skelly001/T1D_Single_Islets.git
cd T1D_Single_Islets

# The renv environment will be automatically set up when you run Snakemake (see below)
```

## Running the Analysis

### Complete Pipeline

The complete analysis pipeline is automated using Snakemake. After setting up the data (see Data Setup above), run:

``` bash
# Run the complete pipeline with automatic environment setup
snakemake --cores 4
```

**Note**: The `Snakefile` automatically: 1. Initializes the `renv` environment with all required packages 2. Executes all analysis stages in the correct order 3. Manages dependencies between analysis steps

### Analysis Stages

The pipeline consists of 10 major stages that recapitulate the analyses presented in the manuscript:

1.  **MSnSet Creation** - Import and format MSstats data into MSnSet objects
2.  **Preprocessing** - Normalization, filtering, batch correction (standard and iBAQ normalization)
3.  **WGCNA** - Weighted Gene Co-expression Network Analysis for each T1D donor to identify protein modules associated with CD3+ infiltration and insulin intensity
4.  **Immune Signature Identification** - Machine learning (random forest with nested cross-validation) to select the top 40-protein IIRS from 329 candidate proteins identified in WGCNA immune-related modules
5.  **Immune Signature Analysis** - UMAP trajectory visualization, snRNA-seq validation, differential expression analysis (LIMMA), pathway enrichment (CAMERA-PR), and identification of key pathways including:
    -   Antigen processing and presentation
    -   Interferon signaling (Type I and II)
    -   ECM dysregulation (lower association in stage 1 T1D)
    -   Hyaluronan metabolic and glycosaminoglycan catabolic processes
    -   IL-10 production
6.  **Beta Cell Profile** - Identification of the 42-protein BCP by selecting proteins with highest correlation with INS and ENTPD3
7.  **Beta Cell Profile Analysis** - Differential expression (LIMMA), pathway enrichment (CAMERA-PR) revealing:
    -   Mitochondrial translation and gene expression (positive correlation with BCP)
    -   mRNA processing and RNA-splicing (negative correlation with BCP in stage 1 T1D)
    -   ECM remodeling (negative correlation with BCP)
    -   IIRS-BCP trajectory comparison (weak correlation, r = 0.13)
8.  **Clustering** - Fine-grained hierarchical clustering identifying 75 protein modules and their functional associations
9.  **Cell Type QC** - Islet cell type marker validation using scRNA-seq reference data (optional, requires Azimuth preprocessing)
10. **Miscellaneous QC** - Final quality control plots (insulin/glucagon distributions, observed proteins) and donor information summary tables

### Cell Type QC Analysis (Optional)

If you wish to run the cell type quality control analysis (Stage 9), you must first process the Azimuth human pancreas reference dataset:

``` bash
# Navigate to the Azimuth reference directory
cd azimuth-references/human_pancreas_snakemake

# Run the Azimuth pancreas Snakefile
snakemake --cores 4

# Return to main directory
cd ../..

# Then run scripts 9a_1, 9b_1, and 9c_1 manually
Rscript 9a_1-cell_type_markers_prep.R
Rscript 9b_1-cell_type_markers_selection.R
Rscript 9c_1-cell_type_QC_barplot.R
```

The Azimuth reference processing creates the scRNA-seq reference object needed for cell type marker validation.

**Reference**: Azimuth pancreas reference sourced from [satijalab/azimuth-references](https://github.com/satijalab/azimuth-references)

### Running Individual Scripts

You can also run individual analysis scripts manually. If running scripts outside of Snakemake, you must first set up the `renv` environment:

``` bash
# First time setup: Initialize renv environment
Rscript -e "renv::restore(prompt = FALSE)"

# Uncomment the source line in .Rprofile to activate renv for manual script execution
# Edit .Rprofile and uncomment: source("renv/activate.R")

# Then run individual scripts
# Example: Run immune signature selection
Rscript 4b_1-immune_signature_selection.R
```

## Output

Analysis outputs are organized in the `output/` directory by analysis stage:

-   `output/RD1-raw_msnsets/` - Raw MSnSet objects
-   `output/RD2-preprocessing/` - Preprocessed data with batch correction
-   `output/RD3-WGCNA/` - WGCNA co-expression networks and cluster ORA
-   `output/RD4-islet_immune_response_signature/` - IIRS proteins and models
-   `output/RD5-islet_immune_response_signature_analysis/` - IIRS analysis results
-   `output/RD6-beta_cell_profile/` - Beta Cell Profile proteins
-   `output/RD7-beta_cell_profile_analysis/` - Beta Cell Profile analysis results
-   `output/RD8-clustering/` - Clustering results
-   `output/RD9-cell_type_marker_QC/` - Cell type QC
-   `output/RD10-misc/` - QC plots and summary tables

Key output files include: - Heatmaps of protein signatures - Volcano plots of differential expression - UMAP trajectory visualizations - Gene set enrichment results (CAMERA-PR analysis) - Quality control plots

## Key Dependencies

Major R packages used in this analysis: - `MSnbase` / `MSnSet.utils` - Proteomics data structures and manipulation - `limma` - Linear modeling for differential expression analysis and CAMERA-PR competitive gene set testing - `WGCNA` - Weighted Gene Co-expression Network Analysis - `ComplexHeatmap` - Advanced heatmap visualization - `tidyverse` - Data manipulation and visualization - `Seurat` - scRNA-seq reference data processing (for cell type QC) - `mlr3verse` - Machine learning ecosystem for IIRS feature selection - `umap` - Dimensionality reduction and trajectory visualization - `clusterProfiler` - Over-representation analysis (ORA)

All package versions are locked and managed by `renv` (specified in `renv.lock`) to ensure reproducibility.

## Citation

If you use this code or data, please cite:

> Kelly SS, Sarkar S, Williams SM, Fu AD, Butterworth EA, Sagendorf TJ, Nierves LA, Kwon Y, Li X, Petyuk VA, Fulcher JM, Chen J, Nakayasu ES, Atkinson MA, Kulkarni RN, Mathews CE, Zhu Y, Campbell-Thompson M, Qian W-J. **Single-Islet Proteomics Maps Pseudo-Temporal Islet Immune Responses and Dysfunction in Stage 1 Type 1 Diabetes.** bioRxiv 2025.11.10.687674; doi: https://doi.org/10.1101/2025.11.10.687674

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

**Corresponding Authors:** - Wei-Jun Qian (weijun.qian\@pnnl.gov) - Pacific Northwest National Laboratory - Martha Campbell-Thompson (thompmc\@pathology.ufl.edu) - University of Florida - Ying Zhu (yingzhupnnl\@gmail.com) - Pacific Northwest National Laboratory\
- Clayton E. Mathews (cxm\@ufl.edu) - University of Florida

## Acknowledgments

We thank the donors and families of the donors for their invaluable contribution to our research and their help to further understand and hopefully cure type 1 diabetes.  This research was performed with the support of the Network for Pancreatic Organ donors with Diabetes (nPOD; RRID:SCR_014641), a collaborative type 1 diabetes research project supported by Breakthrough T1D and The Leona M. & Harry B. Helmsley Charitable Trust (Grant#3-SRA-2023-1417-S-B). The content and views expressed are the responsibility of this article’s authors and do not necessarily reflect the official view of nPOD. Organ Procurement Organizations (OPO) partnering with nPOD to provide research resources are listed at [**https://npod.org/for-partners/npod-partners/**](https://npod.org/for-partners/npod-partners/). This research was supported by NIH Grants R01DK122160, R01DK135081, R01DK131059, R01DK123329, P01AI42288, and U01DK137113. Proteomics was performed in the Environmental Molecular Sciences Laboratory, a national scientific user facility sponsored by the DOE and located at Pacific Northwest National Laboratory, which is operated by Battelle Memorial Institute for the DOE under Contract DE-AC05-76RL0 1830.  This work utilized a LEICA 7000 laser microdissection microscope purchased with a NIH shared instrumentation grant S10OD016350 and operated by the University of Florida Molecular Pathology Core (RRID:SCR_016601).