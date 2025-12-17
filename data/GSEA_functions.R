# ==============================================================================
# GSEA Functions
# ==============================================================================
# Description: Helper functions for Gene Set Enrichment Analysis (GSEA) using
#              clusterProfiler and ReactomePA. Includes functions for creating
#              gene lists, performing GSEA across multiple databases (GO,
#              Reactome, KEGG), simplifying redundant terms, generating
#              visualizations (ridgeplots, heatmaps), and creating KEGG pathway
#              diagrams.
#
# Usage: Source this file to load GSEA utility functions for downstream
#        enrichment analyses.
# ==============================================================================

# Imported by gseGO
library(AnnotationDbi)
library(BiocGenerics)
library(stats4)
library(parallel)
# Imported by gsePathway
library(ReactomePA)
# Others
library(clusterProfiler) # GSEA
library(dplyr)
library(data.table) # rbindlist()
library(org.Hs.eg.db) # Human database
library(ggplot2)
library(scales)
library(rlang)
library(stringr)
library(stats)
library(ggridges)
library(reactome.db)
# For heatmaps
library(tidyselect)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
# KEGG
library(pathview) # KEGG diagrams
library(xml2) # Fix KEGG .xml files


## Create geneList for GSEA ----
# Creates a named vector of ranking metric values for GSEA.
# The values are -log10(p-value) * sign(logFC) and the names are
# the Entrez gene IDs.
create_genelist <- function(df, level,
                            metric = "-log10(P.Value) * sign(logFC)")
{
  if (timewise) {
    df_sub <- filter(df, sex_training == level, !is.na(entrez_id))
  } else {
    df_sub <- filter(df, training_group == level, !is.na(entrez_id))
  }

  # If there are duplicates of the same entrez_id,
  # select the row with the smallest absolute ranking metric value
  geneList <- df_sub %>%
    mutate(ranking_metric = !!parse_expr(metric),
           entrez_id = as.character(entrez_id)) %>%
    group_by(entrez_id) %>%
    # If there are multiple entries for a particular Entrez ID,
    # take the average of the ranking metric.
    summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>%
    arrange(-ranking_metric) %>% # Sort for GSEA
    filter(!is.na(ranking_metric)) %>%
    # Ranking metric vector. Can include ties, but increase nPermSimple
    deframe() # Convert to named vector

  return(geneList)
}


## Generate rds files of GSEA results ----
# We use the signed log p-value as the ranking metric. This is similar
# to using the t-statistic, but the values have a higher peak, so the results
# may be more conservative. A density plot of the signed log p-values
# compared to that of the t-statistics would show this.
create_GSEA_data <- function(df, levels,
                             metric = "-log10(P.Value) * sign(logFC)",
                             timewise = TRUE,
                             database = c("BP", "CC", "MF",
                                          "Reactome", "KEGG"),
                             new_args = list())
{
  for (db in database) {
    message(db)

    for (level in levels) {
      message(paste0("  ", level))

      # Create GSEA input vector
      geneList <- create_genelist(df = df, level = level, metric = metric)

      # GSEA ----
      # List of arguments that are used for all GSEA functions
      common_args <- list(geneList = geneList,
                          # Do not filter by p-values
                          pvalueCutoff = 1,
                          eps = 0,
                          verbose = FALSE,
                          pAdjustMethod = "BH",
                          nPermSimple = 10000)

      # Add additional arguments and determine GSEA function to use
      if (db %in% c("BP", "CC", "MF")) {
        gsea_args <- c(list(ont = db,
                            OrgDb = org.Hs.eg.db,
                            keyType = "ENTREZID"),
                       common_args)
        gsea_fun <- gseGO
      } else if (db == "Reactome") {
        gsea_args <- c(list(organism = "human"), common_args)
        gsea_fun <- gsePathway
      } else if (db == "KEGG") {
        gsea_args <- c(list(organism = "hsa"), common_args)
        gsea_fun <- gseKEGG
      }

      # Update arguments from user
      gsea_args <- modifyList(x = gsea_args,
                              val = new_args,
                              keep.null = TRUE)

      # Perform GSEA
      x <- do.call(what = gsea_fun, args = gsea_args)

      # Remove gene sets that were not annotated to reduce file size
      x@geneSets <- x@geneSets[names(x@geneSets) %in% x@result$ID]

      message(sprintf("    Result has %d rows", nrow(x@result)))

      # Create data folder if it doesn't exist
      if (!dir.exists("data")) {
        dir.create("data")
      }

      # Create data and RDS_files folders if they don't exist
      if (!dir.exists("data/RDS_files")) {
        dir.create("data/RDS_files")
      }

      # Dynamic file name
      file_name <- ifelse(timewise, "GSEA_%s_%s_vs_SED.rds",
                          "GSEA_%s_MvF_%s.rds")
      file_name <- sprintf(file.path("data/RDS_files", file_name), db, level)

      if (nrow(x@result >= 1)) {
        # Save enrichment results
        saveRDS(x, file = file_name)
      }
    }
  }
}



## Create lists to simplify GO terms with simplify_to_general_GO ----
# onts <- c("BP", "CC", "MF")
# # Create a list of 3 lists
# ont_list <- lapply(onts, function(ont) {
#   # For levels 1-14, get the GO terms
#   lapply(1:14, function(level) {
#     clusterProfiler:::getGOLevel(ont, level)
#   })
# })
# # Set the names to easily access elements later
# names(ont_list) <- onts

# Take a look at the results
# lapply(ont_list, lapply, head) # first 6 terms in each level




## Simplify to parent GO terms ----
# If two GO terms have a semantic similarity > 0.7, the term with the lower
# level (the more general a term, the lower the level) is retained.
# If two terms have the same level, both are retained.
simplify_to_general_GO <- function(x, cutoff = 0.7) {

  if (x@setType %in% c("BP", "MF", "CC")) {
    # For each ID, get the levels where it appears.
    # Also get the first level where it appears.
    x@result <- x@result %>%
      mutate(
        # All levels
        GO_levels = sapply(ID, function(ID) {
          which(unlist(
            lapply(ont_list[[x@setType]], function(s) ID %in% s)
          ))
        }),
        # First level
        first_GO_level = unlist(lapply(GO_levels, min)),
        # Convert list to character vector
        GO_levels = unlist(lapply(GO_levels, function(i) {
          paste(unlist(i), collapse = ", ")
        }))
      )

    # If two terms have a similarity > 0.7, select the term in the most
    # general level (lowest level number).
    # This will likely select the parent term.
    x <- simplify(x, cutoff = cutoff, by = "first_GO_level", select_fun = min)

    # Remove unnecessary columns
    x@result <- x@result %>%
      select(-GO_levels, -first_GO_level)
  }

  return(x)
}



## Combine all GSEA .rds files into a single .txt file ----
# We also simplify to parent GO terms.
create_combined_GSEA_table <- function(files, timewise = TRUE,
                                       ome = "global_proteomics") {

  res <- lapply(files, function(file) {

    x <- readRDS(file)

    if (timewise) {
      # Add sex column to timewise results
      x@result$sex <- sex <- ifelse(grepl("female", file), "female", "male")
      training_group <- gsub(".*male_(.*)_vs.*", "\\1", file)
      db <- gsub("data/RDS_files/GSEA_([^_]+).*", "\\1", file)

      message(sprintf("%s_%s %s", sex, training_group, db))
    } else {
      db <- gsub("data/RDS_files/GSEA_(.*)_MvF.*", "\\1", file)
      training_group <- gsub(".*_(.*)\\.rds", "\\1", file)

      message(paste(training_group, db))
    }

    # Filter to significant sets and add new columns
    x@result <- x@result %>%
      filter(p.adjust < 0.05) %>%
      mutate(database = db,
             training_group = training_group)

    # Continue if there are entries after filtering
    if (nrow(x@result) >= 1) {
      ## Handling redundant GO terms
      if (db %in% c("BP", "MF", "CC")) {
        x <- simplify_to_general_GO(x, cutoff = 0.7)
      }

      # Map entrez IDs to gene symbols
      x <- setReadable(x, "org.Hs.eg.db", keyType = "ENTREZID")

      x@result
    }
  }) %>%
    rbindlist(fill = TRUE) # Combine all tables in the list

  # Save full table
  write.table(res, file = "data/combined_GSEA_results.txt",
              quote = FALSE, row.names = FALSE, sep = "\t")
}




## Format GSEA term descriptions for plotting ----
# Capitalize gene set descriptions and reduce length of long descriptions
format_description <- function(x, max_len = 60) {
  # Capitalize first letter of each description
  # (unless it begins with certain exceptions)
  exceptions <- c("[^ -]+[RD]NA", "cGMP")

  x <- x %>%
    mutate(Description = ifelse(
      grepl(paste0("^", paste(exceptions, collapse = "|^")), Description),
      Description,
      paste0(toupper(substr(Description, 1, 1)),
             substr(Description, 2, nchar(Description))))
    )

  # If Description is too long, cut and add ... (set ID)
  for (i in 1:nrow(x)) {
    x$Description[i] <- ifelse(
      nchar(x$Description[i]) > max_len,
      sprintf("%s... (%s)",
              sapply(x$Description[i],
                     function(des) {
                       strwrap(des, width = max_len - (nchar(x$ID[i]) + 10))
                     }), x$ID[i]),
      x$Description[i])
  }

  return(x)
}




## Automatically determine best transformation for plotting ----
# Determine whether the identity or log10 transform is better.
# If the powers in scientific notation are at most 1 away from each other,
# Use the identity transformation. Otherwise, log10.
# Also generates nice breaks for plotting.
auto_transform <- function(values, n = 4) {

  powers <- as.numeric(gsub(".*e([^e])", "\\1", scientific(values)))

  if (diff(range(powers)) <= 2) {
    trans <- "identity"
    breaks <- pretty_breaks(n = n)
  } else {
    trans <- "log10"
    breaks <- trans_breaks("log10", function(x) 10 ^ x, n = n)
  }

  return(list(trans = trans, breaks = breaks))
}




#' @title Create GSEA ridgeplots
#'
#' @description This is a modified version of clusterProfiler::ridgeplot that
#'   makes use of `auto_transform` and pretty breaks.
#'
#' @param x GSEA result object produced by clusterProfiler or ReactomePA.
#' @param showCategory numeric; number of most significantly-enriched terms to plot.
#' @param fill character; name of a column in `x@result` used for the ridgeplot fill.
#' @param core_enrichment logical; whether to only plot the density of the core enrichment genes. Default is `TRUE`.
#' @param orderBy character; name of column in `x@result` used to order the terms. Default is "NES", which orders the y-axis by NES.
#' @param decreasing logical; determines the direction of the `orderBy` column. The default (`decreasing = FALSE`) orders the y-axis in descending order.

ridgeplot_new <- function(x,
                          showCategory = 30,
                          fill = "p.adjust",
                          core_enrichment = TRUE,
                          orderBy = "NES",
                          decreasing = FALSE)
{
  n <- showCategory

  n_terms <- which(unlist(lapply(geneInCategory(x), length)) >= 4)

  if (core_enrichment) {
    gs2id <- geneInCategory(x)[n_terms[seq_len(n)]]
  }
  else {
    gs2id <- x@geneSets[x$ID[seq_len(n)]]
  }
  gs2val <- lapply(gs2id, function(id) {
    res <- x@geneList[id]
    res <- res[!is.na(res)]
  })
  gs2val <- gs2val[lapply(gs2val, length) > 0]
  nn <- names(gs2val)
  i <- match(nn, x$ID)
  nn <- x$Description[i]

  j <- order(x@result[[orderBy]][i], decreasing = decreasing)

  len <- sapply(gs2val, length)
  gs2val.df <- data.frame(category = rep(nn, times = len),
                          color = rep(x[i, fill], times = len),
                          value = unlist(gs2val))

  # Choose transform and breaks
  y <- auto_transform(gs2val.df$fill)

  colnames(gs2val.df)[2] <- fill

  gs2val.df$category <- factor(gs2val.df$category, levels = nn[j])

  ggplot(gs2val.df,
         aes_string(x = "value", y = "category", fill = fill)) +
    geom_density_ridges(size = 0.4) +
    scale_fill_gradient(name = "Adjusted\np-value",
                        trans = y$trans,
                        breaks = y$breaks,
                        low = "darkred", high = "grey90",
                        labels = label_scientific(digits = 2))
}




#' @title Generate GSEA ridgeplots for each .rds file
#'
#' @description Uses `ridgeplot_new` to construct ridgeplots for all comparisons
#'   and databases.
#'
#' @param files character; RDS files containing clusterProfiler GSEA result
#'   objects.
#' @param timewise logical; whether the results were from timewise comparisons.
#' @param showCategory numeric; the number of most significantly-enriched terms
#'   to display in the ridgeplot. Default is 30.
#' @param ome character; ome used as part of the output file name.
#' @param xlab character; title of x-axis.
#' @param file_type character; type of file to save the ridgeplot as. Default is
#'   "png".
#' @param dpi numeric; file resolution (dots per inch). Default is 300
#'   (standard).

create_GSEA_ridgeplots <- function(files, timewise = TRUE, showCategory = 30,
                                   ome = "proteomics",
                                   xlab = "Gene Ranking Metric",
                                   file_type = "png",
                                   dpi = 300)
{
  for (f in files) {
    if (timewise) {
      sex <- ifelse(grepl("female", f), "female", "male")
      training_group <- gsub(".*male_(.*)_vs.*", "\\1", f)
      db <- gsub("data/RDS_files/GSEA_([^_]+).*", "\\1", f)

      db_long <- switch(db,
                        "BP" = "Biological Processes",
                        "CC" = "Cellular Components",
                        "MF" = "Molecular Functions",
                        "Reactome" = "Reactome Pathways")

      # Dynamic plot title
      plot_title <-
        sprintf("Most Significantly-Enriched %s:\n %s wk Relative to Baseline %ss",
                db_long, substring(training_group, 1, 1), str_to_title(sex))

      # Dynamic file name
      file_name <- sprintf("ridgeplots/%s/%s_timewise_GSEA_%s_%s_%s",
                           db, ome, db, sex, training_group)

      message(sprintf("%s_%s %s", sex, training_group, db))

    } else {
      db <- gsub("data/RDS_files/GSEA_(.*)_MvF.*", "\\1", f)
      training_group <- gsub(".*_(.*)\\.rds", "\\1", f)

      db_long <- switch(db,
                        "BP" = "Biological Processes",
                        "CC" = "Cellular Components",
                        "MF" = "Molecular Functions",
                        "Reactome" = "Reactome Pathways")

      # Dynamic plot title

      tg <- ifelse(grepl("[1248]", training_group),
                   paste(substr(training_group, 1, 1),
                         ifelse(grepl("1", training_group),
                                "wk", "wks")),
                   "Baseline")

      plot_title <- sprintf("Most Significantly-Enriched %s:\nMales Relative to Females at %s",
                            db_long, tg)

      # Dynamic file name
      file_name <- sprintf("ridgeplots/%s/%s_MvF_GSEA_%s_%s",
                           db, ome, db, training_group)

      message(paste(training_group, db))
    }

    x <- readRDS(f)
    x@result <- filter(x@result, p.adjust < 0.05)

    # Continue if there are entries after filtering
    if (nrow(x@result) >= 1) {
      ## Handling redundant GO terms
      if (db != "Reactome") {
        x <- simplify_to_general_GO(x, cutoff = 0.7)
      }

      # Format description
      x@result <- format_description(x@result)

      # Ridgeplot
      p <- ridgeplot_new(x, showCategory = showCategory)

      # quant <- 0.999
      quant <- 1
      # Middle 99.9% of data
      xlims <- quantile(p$data$value, (1 - quant) / 2 + c(0, quant))

      p <- p +
        scale_x_continuous(xlab, limits = xlims, expand = expansion(mult = 0)) +
        coord_cartesian(clip = "off", xlim = xlims) +
        labs(title = plot_title,
             subtitle = "Distribution of Core Enrichment Genes") +
        scale_y_discrete(name = NULL) +
        geom_vline(xintercept = 0, lty = "longdash") +
        guides(fill = guide_colorbar(barwidth = 0.5,
                                     barheight = 5)) +
        theme_bw(base_size = 0.7*12) +
        theme(
          axis.ticks = element_line(size = rel(1), color = "black"),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid.major.y = element_line(color = "black"),
          axis.line = element_line(color = "black"),
          text = element_text(color = "black"),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(lineheight = 0.875,
                                     margin = margin(8, 4, 8, 0, "pt")),
          legend.title = element_text(hjust = 0.5, vjust = 0,
                                      margin = margin(0, 0, 4, 0, "pt")),
          plot.title = element_text(size = rel(1.15),
                                    hjust = 0.5,
                                    margin = margin(0)),
          plot.subtitle = element_text(hjust = 0.5,
                                       color = "grey35",
                                       margin = margin(2, 0, 2, 0, "pt")),
          legend.text = element_text(hjust = 0),
          axis.title.x = element_text(margin = margin(4, 0, 0, 0, "pt")),
          legend.margin = margin(0))

      # Create plots and db folders if they don't exist
      if (!dir.exists("ridgeplots")) {
        dir.create("ridgeplots")
      }

      if (!dir.exists(file.path("ridgeplots", db))) {
        dir.create(file.path("ridgeplots", db))
      }

      # Save plot with dynamic file name
      ggsave(paste0(file_name, ".", file_type), p,
             # Dynamic plot height
             height = 2.5 + 0.042 * length(unique(p$data$category)),
             width = 5.2, dpi = dpi)
    }
  }
}




#' @title Create input file for create_GSEA_heatmap
#'
#' @description Reformats the combined GSEA file into a matrix with columns for
#'   the adjusted p-values and NES (1 column each per comparison). For MvF
#'   results, this means 5 columns each for p-values and NES. For timewise, it
#'   is 8 columns each.
#'
#' @param combined_file character; name of file containing all significant GSEA
#'   results for a particular ome/tissue. This is used in conjunction with
#'   `min_count` to select potential terms.
#' @param files character; names of RDS files containing enrichResult objects
#'   produced by clusterProfiler.
#' @param timewise logical; whether the input is from timewise comparisons.
#' @param min_count numeric; minimum number of comparisons where a term shows up
#'   as significantly enriched. Default 1 comparison.
#' @param file_name character; output filename

create_heatmap_data <- function(combined_file, # Combined file
                                files, # Vector of files to combine
                                timewise, # logical
                                min_count = 1,
                                file_name = "data/GSEA_full_heatmap_matrix.txt") {

  # Get the IDs from the combined_file. These are the significant and
  # more broad terms.
  temp <- read.delim(combined_file)
  terms <- temp %>%
    group_by(ID) %>%
    tally() %>%
    filter(n >= min_count) %>%
    pull(ID)

  ## Combine all .rds files and filter to ID in terms
  res <- lapply(files, function(f) {

    x <- readRDS(f)

    if (timewise) {
      sex <- ifelse(grepl("female", f), "female", "male")
      training_group <- gsub(".*male_(.*)_vs.*", "\\1", f)
      db <- gsub("data/RDS_files/GSEA_([^_]+).*", "\\1", f)

      message(sprintf("%s_%s %s", sex, training_group, db))

      # Add sex column to timewise results
      x@result <- x@result %>%
        mutate(sex = factor(sex, levels = c("female", "male")))

    } else {
      db <- gsub("data/RDS_files/GSEA_(.*)_MvF.*", "\\1", f)
      training_group <- gsub(".*_(.*)\\.rds", "\\1", f)

      message(paste(training_group, db))
    }

    x@result %>%
      mutate(database = db,
             training_group = training_group)
  }) %>%
    rbindlist(fill = TRUE) %>%
    filter(ID %in% terms) %>%
    select(ID, Description, NES, p.adjust,
           database, training_group, any_of("sex")) %>%
    mutate(training_group = factor(
      training_group,
      levels = {if ("sex" %in% names(.))
        paste0(2 ^ (0:3), "W")
        else c("SED", paste0(2 ^ (0:3), "W"))})) %>%
    arrange({if ("sex" %in% names(.)) sex else NULL},
            training_group)

  # Table of adjusted p-values
  p_vals <- res %>%
    select(-NES) %>%
    pivot_wider(values_from = p.adjust,
                names_from = c(training_group, any_of("sex")),
                names_prefix = "pvals_")

  # Table of normalized enrichment scores
  NES <- res %>%
    select(-p.adjust) %>%
    pivot_wider(values_from = NES,
                names_from = c(training_group, any_of("sex")),
                names_prefix = "NES_")

  # Combine p_vals and NES
  res <- left_join(p_vals, NES)

  # Save table
  write.table(res, file = file_name, sep = "\t",
              row.names = FALSE, quote = FALSE)
}




#' @title Create heatmap for each database
#'
#' @description Leverages the ComplexHeatmap package to create a dot heatmap
#'   where the color of the circles is determined by the NES, the size is
#'   determined by -log10(adjusted p-value) (larger dots mean more significant),
#'   and the background fill is based on whether the adjusted p-value passes the
#'   0.05 threshold. The rows are clustered using the Pearson correlation of the
#'   NES.
#'
#' @param file character; name of file produced by `create_heatmap_data`.
#' @param n_top numeric; number of most significant (on-average) terms to
#'   display in the heatmap. Default is 40.
#' @param tissue character; tissue included in the heatmap file name.
#' @param ome character; ome included in the heatmap file name.
#' @param database character; one or more of "BP", "MF", "CC", "Reactome", and
#'   "KEGG". All databases are considered, by default.
#' @param terms character; term IDs. If provided, these will be displayed in the
#'   heatmap instead of the `n_top` terms.
#' @param suffix character; suffix for file name. Default is an empty string "".
#' @param scale_by_row logical; whether to scale the sizes of the circles by row
#'   such that the most significant comparison is the point of maximum size.
#'   This removes the ability to compare the sizes of circles between rows
#'   (though we can still compare the colors), but it is easier to see term-wise
#'   trajectories.
#' @param z numeric; scaling factor for plot. All heatmap elements will be
#'   scaled proportionately. All this really does is increase the resolution and
#'   size of the final heatmap.

create_GSEA_heatmaps <- function(file, n_top = 40, tissue, ome,
                                 database = c("BP", "MF", "CC",
                                              "Reactome", "KEGG"),
                                 terms = character(0),
                                 suffix = "", scale_by_row = TRUE, z = 1)
{
  # Input file produced by create_heatmap_data
  input <- read.delim(file)

  # Loop over each specified database
  for (db in database) {
    # Filter to specific database and format Description
    input_db <- filter(input, database == db) %>%
      format_description(max_len = 40)

    # Filter to specific terms, if provided
    if (!identical(terms, character(0))) {
      input_db <- filter(input_db, ID %in% terms)
      n_top <- nrow(input_db)
    }

    # Matrix of adjusted p-values only
    pval_mat <- input_db %>%
      select(Description, starts_with("pvals_")) %>%
      column_to_rownames("Description") %>%
      as.matrix()

    # If there are 8 columns, the results are timewise
    timewise <- ncol(pval_mat) == 8

    if (timewise) {
      timepoint_rows <- 2
      by_row <- TRUE
    } else {
      timepoint_rows <- 1
      by_row <- FALSE
    }

    # Average adjusted p-value for each term
    mean_pval <- apply(pval_mat, 1, mean)

    # Select the n_top terms with the lowest average adjusted p-value
    pval_mat <- pval_mat[order(mean_pval, decreasing = FALSE), ] %>%
      .[1:min(n_top, nrow(pval_mat)), ]

    # Matrix of normalized enrichment scores only
    NES_mat <- input_db %>%
      select(Description, starts_with("NES")) %>%
      filter(Description %in% rownames(pval_mat)) %>%
      column_to_rownames("Description") %>%
      as.matrix()

    # Make sure row order is the same as pval_mat
    NES_mat <- NES_mat[rownames(pval_mat), ]

    # Save matrices
    write.table(rownames_to_column(
      as.data.frame(cbind(pval_mat, NES_mat)), "feature"
    ),
    file = sprintf("data/GSEA_heatmap_matrix_%s.txt", db),
    quote = FALSE, row.names = FALSE, sep = "\t")

    # Color function
    if (all(sign(NES_mat) %in% c(0, +1))) {
      col_breaks <- c(0, max(NES_mat, na.rm = TRUE))
      col_colors <- c("white", "red")
      temp_fun <- circlize::colorRamp2(breaks = col_breaks,
                                       colors = col_colors)
      col_breaks <- c(0, max(NES_mat, na.rm = TRUE) * c(1, 0.5))
      col_colors <- temp_fun(col_breaks)
    } else if (all(sign(NES_mat) %in% c(0, -1))) {
      col_breaks <- c(min(NES_mat, na.rm = TRUE), 0)
      col_colors <- c("blue", "white")
      temp_fun <- circlize::colorRamp2(breaks = col_breaks,
                                       colors = col_colors)
      col_breaks <- c(min(NES_mat, na.rm = TRUE) * c(0.5, 1), 0)
      col_colors <- temp_fun(col_breaks)
      col_colors <- sub("FF$", "", col_colors)
    } else {
      col_breaks <- c(min(NES_mat), 0, max(NES_mat))
      col_colors <- c("blue", "white", "red")
    }

    col_fun <- circlize::colorRamp2(breaks = col_breaks,
                                    colors = col_colors)

    # Colors for column annotation
    colors <- list(Training = c("#F7FCB9", "#ADDD8E", "#238443", "#002612"))
    names(colors[["Training"]]) <- c("1W", "2W", "4W", "8W")

    # If not timewise (MvF), add SED color
    if (!timewise) {
      colors[["Training"]] <- c("white", colors[["Training"]])
      names(colors[["Training"]])[1] <- "SED"
    }

    # Column annotation data frame
    anno_df <- data.frame(Training = rep(names(colors[["Training"]]),
                                         ifelse(timewise, 2, 1))) %>%
      mutate(Training = factor(Training, levels = unique(Training)))

    Training_list <- list(direction = "horizontal",
                          nrow = 2,
                          by_row = T,
                          grid_height = unit(6, "points") * z,
                          grid_width = unit(6, "points") * z,
                          border = "black",
                          labels_gp = gpar(fontsize = 6 * z),
                          title_gp = gpar(fontsize = 6 * z))

    # If timewise, add Sex list component and data frame column
    if (timewise) {
      colors[["Sex"]] <- c("#ff6eff", "#5555ff")
      names(colors$Sex) <- c("Female", "Male")

      anno_df <- anno_df %>%
        mutate(Sex = rep(c("Female", "Male"), each = 4),
               Sex = factor(Sex, levels = unique(Sex)))

      # Column annotation
      ha <- HeatmapAnnotation(
        df = anno_df, which = "column", gp = gpar(col = "black"),
        col = colors, gap = 0, show_annotation_name = FALSE,
        annotation_legend_param = list(
          Sex = list(direction = "horizontal",
                     # nrow = 1,
                     ncol = 1,
                     nrow = timepoint_rows,
                     grid_height = unit(6, "points") * z,
                     grid_width = unit(6, "points") * z,
                     border = "black", # Border around each box
                     labels_gp = gpar(fontsize = 6 * z),
                     title_gp = gpar(fontsize = 6 * z)
          ),
          Training = Training_list
        )
      )
    } else {
      # Column annotation
      ha <- HeatmapAnnotation(
        df = anno_df, which = "column", gp = gpar(col = "black"),
        col = colors, gap = 0, show_annotation_name = FALSE,
        annotation_legend_param = list(
          Training = Training_list
        )
      )
    }

    # Row annotation
    ra <- rowAnnotation(terms = anno_text(
      x = rownames(NES_mat),
      gp = gpar(col = "black",
                fontsize = 6 * z
      ))
    )

    # Heatmap dimensions. Each cell is a square with width = height = a
    a <- unit(7, "points") * z
    width <- ncol(pval_mat) * a
    height <- nrow(pval_mat) * a

    db_long <- case_when(db == "BP" ~ "Biological Processes",
                         db == "CC" ~ "Cellular Components",
                         db == "MF" ~ "Molecular Functions",
                         TRUE ~ "Reactome Pathways")
    # column_title <- paste("Top Enriched", db_long)
    ome_new <- case_when(ome == "RNA-Seq" ~ "TRNSCRPT",
                         ome == "proteomics" ~ "PROT",
                         TRUE ~ "PHOSPHO")
    db_new <- ifelse(db %in% c("BP", "MF", "CC"),
                     paste0("GO-", db),
                     db)
    column_title <- paste(ome_new, "GSEA:", db_new)

    # Create heatmap.
    # Each cell is defined by cell_fun, which places a grey border around
    # each cell and adds a filled circle with fill determined by NES and
    # size determined by -log10(unadjusted p-value). The size of circles is
    # scaled by row so that the most significant comparison group for
    # each gene set has maximum size. This is good when there are terms that
    # are extremely significant (outliers), as they would make all other
    # circles difficult to see; however, it we are not able to directly
    # compare rows, and there is no way to include an overall size legend.
    # Since the rows of the heatmap are clustered using the Pearson correlation
    # of the NES, we are still able to see which terms are similar.
    ht <- Heatmap(NES_mat, # Dendrogram based on NES
                  col = col_fun,
                  border = "black",
                  heatmap_legend_param = list(
                    at = round(col_breaks, 1),
                    title = "NES",
                    border = "black",
                    title_gp = gpar(fontsize = 6 * z),
                    labels_gp = gpar(fontsize = 6 * z),
                    legend_direction = "horizontal",
                    legend_width = unit(0.4, "in") * z
                  ),
                  width = width,
                  height = height,
                  top_annotation = ha,
                  # column_title = column_title, #
                  column_title_gp = gpar(fontsize = 8 * z),
                  cluster_columns = FALSE,
                  show_column_names = FALSE,
                  rect_gp = gpar(fill = NA, col = NA),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width, height,
                              gp = gpar(col = NA,
                                        fill = ifelse(pval_mat[i, j] < 0.05,
                                                      "lightgrey", NA)))
                    if (scale_by_row) {
                      grid.circle(
                        x = x,
                        y = y,
                        # Radius scaled by row-wise max value
                        # a / 2 is half of the cell dimensions
                        # Radius is at least 0.15 of max cell dim
                        r = max(-log10(pval_mat[i, j])/
                                  max(-log10(pval_mat[i, ])), 0.15)*a/2,
                        gp = gpar(fill = col_fun(NES_mat[i, j]),
                                  col = NA)
                      )
                    } else {
                      grid.circle(x = x,
                                  y = y,
                                  # Radius scaled by max value of matrix
                                  r = -log10(pval_mat[i, j]) /
                                    max(-log10(pval_mat)) * a / 2,
                                  gp = gpar(fill = col_fun(NES_mat[i, j]),
                                            col = NA))
                    }
                  }) +
      ra

    # Legend for background fill
    lt <- Legend(at = 1:2,
                 labels = c("< 0.05", expression("">="0.05")),
                 title = "BH-adjusted\np-value",
                 labels_gp = gpar(fontsize = 6 * z),
                 title_gp = gpar(fontsize = 6 * z),
                 legend_gp = gpar(fill = c("lightgrey", "white")),
                 border = "black", direction = "horizontal", nrow = 2,
                 grid_height = unit(6, "points") * z,
                 grid_width = unit(6, "points") * z
    )

    # Create heatmaps folder
    if(!dir.exists("heatmaps")) {
      dir.create("heatmaps")
    }

    # Save heatmap
    tiff(sprintf("heatmaps/heatmap_%s_%s_%s_%s%s.tiff",
                 tissue, ome, ifelse(timewise, "timewise", "MvF"), db, suffix),
         compression = "lzw",
         height = 1.4 + 4.1 / 40 * n_top * z,
         width = 2.8 * z,
         units = "in",
         res = 300)
    draw(ht,
         heatmap_legend_list = list(lt),
         merge_legends = TRUE,
         heatmap_legend_side = "bottom",
         legend_gap = unit(0.16, "in") * z # space between legends
    )
    dev.off()
  }
}



#' @title Create KEGG Pathway Diagrams
#'
#' @description First, downloads and modifies KEGG files for all significant
#'   pathways. Then, uses `pathview` to construct diagrams--saving them as .png
#'   files in the current working directory. Finally, files are moved to a
#'   "kegg_pathview_diagrams" folder.
#'
#' @param df
#' @param timewise logical;
#' @param metric character; the ranking metric used to determine the range of
#'   values for the pathview color scale.
#' @param download_xml logical; whether to download the KEGG xml files.

kegg_pathview <- function(df, timewise = TRUE,
                          metric = "-log10(P.Value) * sign(logFC)",
                          download_xml = TRUE)
{
  files <- list.files("./data/RDS_files", pattern = "KEGG", full.names = TRUE)

  for (file in files) {
    x <- readRDS(file)

    if (timewise) {
      # Add sex column to timewise results
      x@result$sex <- sex <- ifelse(grepl("female", file), "female", "male")
      training_group <- gsub(".*male_(.*)_vs.*", "\\1", file)
      sex_training <- paste0(sex, "_", training_group)

      lvl <- sex_training

      message(sprintf("%s_%s KEGG", sex, training_group))
    } else {
      training_group <- gsub(".*_(.*)\\.rds", "\\1", file)
      sex <- NULL
      lvl <- training_group
      message(training_group)
    }

    # Significant KEGG terms
    significant_terms <- x@result$ID[x@result$p.adjust < 0.05]

    if (length(significant_terms) > 0) {
      geneList <- create_genelist(df = df,
                                  metric = metric,
                                  level = lvl)
      limit <- list(gene = round(max(abs(geneList), na.rm = T)), cpd = 1)
      # limit <- list(gene = range(geneList), cpd = 1)

      # File suffix
      out.suffix <- gsub(".*GSEA_KEGG_(.*)\\.rds", "\\1", file)

      if (!dir.exists("./kegg_files"))
        dir.create("./kegg_files")

      if (download_xml) {
        message("Downloading XML and PNG files.")
        kegg.files <- list.files("./kegg_files", pattern = "xml")
        kegg.files <- sub("\\.xml", "", kegg.files)
        pathway.id <- setdiff(significant_terms, kegg.files)

        if (length(pathway.id) >= 1) {
          pathway.id <- sub("hsa", "", significant_terms)
          download.kegg(pathway.id = pathway.id,
                        species = "hsa",
                        kegg.dir = "./kegg_files",
                        file.type = c("xml", "png"))
        }

        message("Removing bad nodes from XML files")
        for (kegg_term in significant_terms) {
          # Fix XML file
          xml_file <- sprintf("./kegg_files/%s.xml", kegg_term)
          x <- read_xml(xml_file)
          all_nodes <- xml_children(x)
          bad_nodes <- xml_find_all(all_nodes, "graphics[@x='1']")
          xml_remove(bad_nodes, free = TRUE) # remove bad nodes
          xml2::write_xml(x, xml_file) # overwrite file
        }
      }

      message("Creating Pathview diagrams")
      ## NOTE: Try the pathview function. If it fails, do not output error.
      ##  Instead, move on to the next file.
      # Pathway diagrams
      # pb <- txtProgressBar(min = 0,
      #                      max = length(significant_terms),
      #                      initial = 0,
      #                      style = 3)
      for (i in seq_along(significant_terms)) {
        try(suppressMessages(
          pathview(gene.data   = geneList,
                   pathway.id  = significant_terms[i],
                   kegg.dir    = "./kegg_files",
                   species     = "hsa",
                   limit       = limit,
                   out.suffix  = out.suffix,
                   na.col      = "grey",
                   low         = list(gene = "blue", cpd = "blue"),
                   mid         = list(gene = "white", cpd = "white"),
                   high        = list(gene = "red", cpd = "red"),
                   map.null    = FALSE)
        ))
      }

      # close(pb) # close progress bar
    }
  }

  ## Remove useless images (no color)
  list.files(pattern = "png", full.names = T) %>%
    {.[!grepl("_", .)]} %>%
    file.remove()

  # Move Pathview diagrams to separate folder
  old_files <- list.files(pattern = "png", full.names = TRUE)
  if (!dir.exists("./kegg_pathview_diagrams"))
    dir.create("./kegg_pathview_diagrams")
  new_files <- sub("\\./", "\\./kegg_pathview_diagrams/", old_files)
  file.rename(from = old_files, to = new_files)
}



