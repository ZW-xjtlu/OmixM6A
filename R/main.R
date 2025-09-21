#' OmixM6A: Methylation State Analysis of m6A Modification Sites
#'
#' This function performs statistical analysis on m6A modification sites using different methods including mixtures of two binomial distributions, beta-binomial distributions, or a mixture of binomial and beta(1,1) distributions. It supports analysis on both individual vectors of counts and `SummarizedExperiment` objects.
#'
#' @param m6A A numeric vector of counts for m6A modifications.
#' @param Total A numeric vector of total counts. Must be the same length as `m6A`.
#' @param se An optional `SummarizedExperiment` object. If provided, the analysis will store results in its assays and metadata.
#' @param method The method for analysis: "bbmix" for beta-binomial mixture, "bmix" for binomial mixture, "zoibmix" for zero-and-one over-inflated binomial mixture, "zoibbmix" for zero-and-one over inflated beta-binomial mixture, "bumix" for binomial and beta(1,1) mixture, "binomial" for simple binomial test. Defaults to "bbmix".
#' @param bbmix_size An integer specifying the maximum number of sites for fitting the model. This parameter is for efficiency considerations and only affects the mixture methods involving beta-binomial. The default value is NULL (no subset).
#' @param cov_threshold An integer specifying a threshold for subsetting sites. Only sites with a total read count greater than or equal to this value are used for model fitting. The default value is 0.
#'
#' @return If `se` is provided, returns a modified `SummarizedExperiment` object with added assays and metadata containing analysis results. Otherwise, returns a data frame with analysis results including m6A counts, total counts, beta values (= m6A/Total), posterior probabilities of foreground (prob_fg), p-values, and BH adjusted p-values (fdr).
#'
#' @details
#' The function calculates several key statistics:
#'
#' - `prob_fg`: The posterior probability of a site being a foreground m6A modification. Useful for predicting m6A sites using a Bayesian classifier approach.
#'
#' - `pvalue`: The posterior predictive p-value, indicating the probability of observing data as or more extreme than the current data under the fitted background model.
#'
#' - `fdr`: The adjusted p-value using the Benjamini-Hochberg (BH) approach, which can control the false discovery rate (fdr).
#'
#' In case of a `SummarizedExperiment` input, the fitted model parameters are stored in its metadata.
#'
#' @import extraDistr SummarizedExperiment S4Vectors MASS
#' @export
#' @examples
#' # Example usage with count vectors
#' library(SummarizedExperiment)
#' m6A_se <- readRDS(system.file("extdata", "example_se.rds", package="OmixM6A"))
#' m6A_counts <- assays(m6A_se)$m6A[,1]
#' total_counts <- assays(m6A_se)$Total[,1]
#' result_df <- OmixM6A(m6A_counts, total_counts)
#' result_df
#'
#' # Example usage with a SummarizedExperiment object
#' result_se <- OmixM6A(se = m6A_se)
#' result_se
#' metadata(result_se) #Check fitted model parameters
#'
OmixM6A <- function(m6A, Total, se = NULL, method = c("bbmix", "bmix", "zoibbmix", "zoibmix", "bumix", "binomial"), bbmix_size = NULL, cov_threshold = 0){
  method <- match.arg(method)
  if(!is.null(se)){
  stopifnot(is(se, "SummarizedExperiment"))
  stopifnot(all(assays(se)$Total-assays(se)$m6A >=0))
  assays(se, withDimnames=FALSE)$beta <- matrix(nrow = nrow(se), ncol = ncol(se))
  assays(se, withDimnames=FALSE)$prob_fg <- matrix(nrow = nrow(se), ncol = ncol(se))
  assays(se, withDimnames=FALSE)$pvalue <- matrix(nrow = nrow(se), ncol = ncol(se))
  assays(se, withDimnames=FALSE)$fdr <- matrix(nrow = nrow(se), ncol = ncol(se))
  metadata(se)$OmixM6A_para <- data.frame(bg_proportion = rep(NA, ncol(se)))
  rownames(metadata(se)$OmixM6A_para) <- colnames(se)
  for(i in seq_len(ncol(se))){
    if(method == "bmix"){
      fit_i <- fit_bmix(assays(se)$m6A[,i], assays(se)$Total[,i], cov_threshold)
    }else if(method == "bbmix"){
      fit_i <- fit_bbmix(assays(se)$m6A[,i], assays(se)$Total[,i], subsize = bbmix_size, cov_threshold)
    }else if(method == "zoibbmix"){
      fit_i <- fit_zoibbmix(assays(se)$m6A[,i], assays(se)$Total[,i], subsize = bbmix_size, cov_threshold)
    }else if(method == "zoibmix"){
      fit_i <- fit_zoibmix(assays(se)$m6A[,i], assays(se)$Total[,i], cov_threshold)
    }else if(method == "bumix"){
      fit_i <- fit_bumix(assays(se)$m6A[,i], assays(se)$Total[,i], cov_threshold)
    }else if(method == "binomial"){
      fit_i <- fit_binomial(assays(se)$m6A[,i], assays(se)$Total[,i], cov_threshold)
    }
    assays(se)$beta[,i] <- fit_i$beta
    assays(se)$prob_fg[,i] <- fit_i$prob_fg
    assays(se)$pvalue[,i] <- fit_i$pvalue
    assays(se)$fdr[,i] <- p.adjust(fit_i$pvalue, method = "fdr")
    for(j in names(fit_i$para)){
      metadata(se)$OmixM6A_para[[j]][i] <- fit_i$para[[j]]
    }
  }
  return(se)
  }else{
    if (method == "bmix") {
      fit <- fit_bmix(m6A, Total, cov_threshold)
    } else if (method == "bbmix") {
      fit <- fit_bbmix(m6A, Total, subsize = bbmix_size, cov_threshold)
    } else if (method == "zoibmix") {
      fit <- fit_zoibmix(m6A, Total, cov_threshold)
    } else if (method == "zoibbmix") {
      fit <- fit_zoibbmix(m6A, Total, subsize = bbmix_size, cov_threshold)
    } else if (method == "bumix"){
      fit <- fit_bumix(m6A, Total, cov_threshold)
    }else if(method == "binomial"){
      fit <- fit_binomial(m6A, Total, cov_threshold)
    }
   result_tbl <- data.frame(m6A_count = m6A,
                            Total_count = Total,
                            beta = fit$beta,
                            prob_fg = fit$prob_fg,
                            pvalue = fit$pvalue,
                            fdr = p.adjust(fit$pvalue, method = "fdr"))
   return(result_tbl)
  }
}

