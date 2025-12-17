# ==============================================================================
# Single Islets Project Functions
# ==============================================================================
# Description: Custom utility functions for processing single islet proteomics
#              data. Includes functions for completeness filtering across donors,
#              imputation using SVD, and renaming tissue block identifiers for
#              visualization.
#
# Usage: Source this file to load project-specific helper functions.
# ==============================================================================

library(MSnSet.utils)
library(dplyr)



# Require xx completeness from every donor --------------------------------

completeness_filt <- \(msnset, threshold) {
    
    compl_thresh <- threshold
    msnset$donor <- droplevels(msnset$donor)
    m_donor <- split(msnset, "donor")
    donor_completeness <- purrr::map(m_donor, \(x){apply(exprs(x), 1, \(xx){mean(!is.na(xx)) >= compl_thresh})})
    
    feat_to_keep <- purrr::map(donor_completeness, tibble::enframe) %>%
        purrr::reduce(\(x, y){left_join(x, y, by = "name")}) %>%
        rowwise() %>%
        mutate(keep_feat = sum(c_across(matches("value"))) == length(m_donor))
    
    m2 <- msnset[feat_to_keep$keep_feat,  ]

}




# Impute ------------------------------------------------------------------

impute_msnset <- \(msnset){
    
    # record feat SD
    stds <- apply(exprs(msnset), 1, \(x){sd(x, na.rm = T)})
    
    # standardize (feat already median centered)
    exprs(msnset) <- t(scale(t(exprs(msnset)), center = F, scale = T))
    
    # split case/control
    m_case <- msnset[, pData(msnset)$class == "case"]
    m_cont <- msnset[, pData(msnset)$class == "control"]
    
    # check for class-dependent missingness
    m_case <- filter_by_occurrence(m_case, 0.001)
    m_cont <- filter_by_occurrence(m_cont, 0.001)
    feat_match <- length(featureNames(m_case)) == length(featureNames(m_cont))
    if (!feat_match) stop("Impute would introduce class-dependent missinness.")
    
    # Impute
    exprs(m_case) <- t(pcaMethods::completeObs(pcaMethods::pca(as(m_case,"ExpressionSet"), method="svdImpute",
                                       nPcs=min(dim(m_case)), center = F)))
    exprs(m_cont) <- t(pcaMethods::completeObs(pcaMethods::pca(as(m_cont,"ExpressionSet"), method="svdImpute",
                                       nPcs=min(dim(m_cont)), center = F)))
    m_imp <- MSnbase::combine(m_case, m_cont)
    
    # ensure feature order matches original
    m_imp <- m_imp[featureNames(msnset), ]
    
    # restore feature magnitude
    exprs(m_imp) <- sweep(exprs(m_imp), 1, STATS = stds, FUN = "*")
    
    return(m_imp)
}




# Rename levels of Block for plotting ------------------------------------

rename_block_for_plot <- \(x){
    m_plot <- msnsets[[x]]
    
    rename_key <- unique(pData(m_plot)$Block) %>% 
        as.character() %>% 
        (\(x) data.frame(old_name = x, new_name = str_extract(x, "[1-9]")))() %>% 
        group_by(new_name) %>% 
        mutate(letter = LETTERS[row_number()],
               new_name2 = paste0(new_name, letter)) %>% 
        pull(new_name2, old_name)
    
    pData(m_plot) <- pData(m_plot) %>% 
        rownames_to_column() %>% 
        mutate(Block = factor(rename_key[Block], levels = sort(rename_key))) %>% 
        column_to_rownames()
    
    m_plot
}


























