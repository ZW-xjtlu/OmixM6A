extract_bbmixFit <- function(m6A, Total, total_threshold = 30, bbmix_size = NULL){
  #require(SummarizedExperiment)
  se <- SummarizedExperiment(assays = list(m6A = cbind(m6A), Total = cbind(Total)), checkDimnames = FALSE)
  indx <- as.vector(assays(se)$Total >= total_threshold)
  N <- sum(indx)
  plot_df <- data.frame(m6A_ratio = rep(NA, 2*N), type = rep(NA, 2*N))
  plot_df$m6A_ratio[1:N] <- assays(se)$m6A[indx]/assays(se)$Total[indx]
  plot_df$type[1:N] <- "estimated from data"
  se_bb <- OmixM6A(se = se[indx,], method = "bbmix", bbmix_size = bbmix_size)
  #Extract the most likely # of bg sites according to the proportion parameter
  prop_fg <- metadata(se_bb)$OmixM6A_para$fg_proportion
  prop_bg <- 1 - prop_fg
  indx_fg_tmp <- order(assays(se_bb)$prob_fg, decreasing = TRUE)[1:round(N*prop_fg)]
  #Convert the index back to dummy variable
  indx_fg <- rep(FALSE, N)
  indx_fg[indx_fg_tmp] <- TRUE
  rm(indx_fg_tmp)
  alpha_bg <- metadata(se_bb)$OmixM6A_para$alpha_m6A_bg
  beta_bg <- metadata(se_bb)$OmixM6A_para$beta_m6A_bg
  alpha_fg <- metadata(se_bb)$OmixM6A_para$alpha_m6A_fg
  beta_fg <- metadata(se_bb)$OmixM6A_para$beta_m6A_fg
  bg_p <- rbeta(sum(indx_fg), alpha_bg, beta_bg)
  fg_p <- rbeta(sum(indx_fg), alpha_fg, beta_fg)
  bg_beta_vals <- rbinom(n = sum(!indx_fg),prob = bg_p, size = assays(se_bb)$Total[!indx_fg])/assays(se_bb)$Total[!indx_fg]
  fg_beta_vals <- rbinom(n = sum(indx_fg),prob = fg_p, size = assays(se_bb)$Total[indx_fg])/assays(se_bb)$Total[indx_fg]
  beta_vals_all <- rep(NA, N)
  beta_vals_all[indx_fg] <- fg_beta_vals
  beta_vals_all[!indx_fg] <- bg_beta_vals
  plot_df$m6A_ratio[(N+1):(2*N)] <- beta_vals_all
  plot_df$type[(N+1):(2*N)] <- "simulated by fitted model"
  qs <- 1-rev(quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], probs = seq(0,1,length.out = 250)))
  qtl_emp <- quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], qs)
  qtl_fit <- quantile(plot_df$m6A_ratio[plot_df$type != "estimated from data"], qs)
  D <- (2-1)+2*2
  BIC <- D*log(N) - 2*sum(log(prop_bg*dbbinom(m6A[indx], Total[indx], alpha_bg, beta_bg) + prop_fg*dbbinom(m6A[indx], Total[indx], alpha_fg, beta_fg)))
  return( list(plot_df_qtl = data.frame(qtl_emp = qtl_emp,
                                          qtl_fit = qtl_fit),
               plot_df_hist = plot_df,
               BIC = BIC) )
}

extract_bumixFit <- function(m6A, Total, total_threshold = 30){
  #require(SummarizedExperiment)
  se <- SummarizedExperiment(assays = list(m6A = cbind(m6A), Total = cbind(Total)), checkDimnames = FALSE)
  indx <- as.vector(assays(se)$Total >= total_threshold)
  N <- sum(indx)
  plot_df <- data.frame(m6A_ratio = rep(NA, 2*N), type = rep(NA, 2*N))
  plot_df$m6A_ratio[1:N] <- assays(se)$m6A[indx]/assays(se)$Total[indx]
  plot_df$type[1:N] <- "estimated from data"
  se_bu <- OmixM6A(se = se[indx,], method = "bumix")
  #Extract the most likely # of bg sites according to the proportion parameter
  prop_fg <- metadata(se_bu)$OmixM6A_para$fg_proportion
  prop_bg <- 1 - prop_fg
  indx_fg_tmp <- order(assays(se_bu)$prob_fg, decreasing = TRUE)[1:round(N*prop_fg)]
  #Convert the index back to dummy variable
  indx_fg <- rep(FALSE, N)
  indx_fg[indx_fg_tmp] <- TRUE
  rm(indx_fg_tmp)
  fg_p <- runif(sum(indx_fg))
  p_bg <- metadata(se_bu)$OmixM6A_para$p_m6A_bg
  bg_beta_vals <- rbinom(n = sum(!indx_fg),prob = p_bg, size = assays(se_bu)$Total[!indx_fg])/assays(se_bu)$Total[!indx_fg]
  fg_beta_vals <- rbinom(n = sum(indx_fg),prob = fg_p, size = assays(se_bu)$Total[indx_fg])/assays(se_bu)$Total[indx_fg]
  beta_vals_all <- rep(NA, N)
  beta_vals_all[indx_fg] <- fg_beta_vals
  beta_vals_all[!indx_fg] <- bg_beta_vals
  plot_df$m6A_ratio[(N+1):(2*N)] <- beta_vals_all
  plot_df$type[(N+1):(2*N)] <- "simulated by fitted model"
  qs <- 1-rev(quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], probs = seq(0,1,length.out = 250)))
  qtl_emp <- quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], qs)
  qtl_fit <- quantile(plot_df$m6A_ratio[plot_df$type != "estimated from data"], qs)
  D <- (2-1)+1
  BIC <- D*log(N) - 2*sum(log(prop_bg*dbinom(m6A[indx], Total[indx], p_bg) + prop_fg*dbbinom(m6A[indx], Total[indx], 1, 1)))
  return( list(plot_df_qtl = data.frame(qtl_emp = qtl_emp,
                                         qtl_fit = qtl_fit),
               plot_df_hist = plot_df,
               BIC = BIC) )
}

extract_bmixFit <- function(m6A, Total, total_threshold = 30){
  #require(SummarizedExperiment)
  se <- SummarizedExperiment(assays = list(m6A = cbind(m6A), Total = cbind(Total)), checkDimnames = FALSE)
  indx <- as.vector(assays(se)$Total >= total_threshold)
  N <- sum(indx)
  plot_df <- data.frame(m6A_ratio = rep(NA, 2*N), type = rep(NA, 2*N))
  plot_df$m6A_ratio[1:N] <- assays(se)$m6A[indx,]/assays(se)$Total[indx,]
  plot_df$type[1:N] <- "estimated from data"
  se_b <- OmixM6A(se = se[indx,], method = "bmix")
  #Extract the most likely # of bg sites according to the proportion parameter
  prop_fg <- metadata(se_b)$OmixM6A_para$fg_proportion
  prop_bg <- 1 - prop_fg
  indx_fg_tmp <- order(assays(se_b)$prob_fg, decreasing = TRUE)[1:round(N*prop_fg)]
  #Convert the index back to dummy variable
  indx_fg <- rep(FALSE, N)
  indx_fg[indx_fg_tmp] <- TRUE
  rm(indx_fg_tmp)
  p_bg <- metadata(se_b)$OmixM6A_para$p_m6A_bg
  p_fg <- metadata(se_b)$OmixM6A_para$p_m6A_fg
  bg_beta_vals <- rbinom(n = sum(!indx_fg),prob = p_bg, size = assays(se_b)$Total[!indx_fg])/assays(se_b)$Total[!indx_fg]
  fg_beta_vals <- rbinom(n = sum(indx_fg),prob = p_fg, size = assays(se_b)$Total[indx_fg])/assays(se_b)$Total[indx_fg]
  beta_vals_all <- rep(NA, N)
  beta_vals_all[indx_fg] <- fg_beta_vals
  beta_vals_all[!indx_fg] <- bg_beta_vals
  plot_df$m6A_ratio[(N+1):(2*N)] <- beta_vals_all
  plot_df$type[(N+1):(2*N)] <- "simulated by fitted model"
  qs <- 1-rev(quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], probs = seq(0,1,length.out = 250)))
  qtl_emp <- quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], qs)
  qtl_fit <- quantile(plot_df$m6A_ratio[plot_df$type != "estimated from data"], qs)
  D <- (2-1)+2
  BIC <- D*log(N) - 2*sum(log(prop_bg*pmax(dbinom(m6A[indx], Total[indx], p_bg), 1e-100) + prop_fg*pmax(dbinom(m6A[indx], Total[indx], p_fg), 1e-100)))
  return( list(plot_df_qtl = data.frame(qtl_emp = qtl_emp,
                                         qtl_fit = qtl_fit),
               plot_df_hist = plot_df,
               BIC = BIC) )
}

extract_binomialFit <- function(m6A, Total, total_threshold = 30){
  #require(SummarizedExperiment)
  se <- SummarizedExperiment(assays = list(m6A = cbind(m6A), Total = cbind(Total)), checkDimnames = FALSE)
  indx <- as.vector(assays(se)$Total >= total_threshold)
  N <- sum(indx)
  plot_df <- data.frame(m6A_ratio = rep(NA, 2*N), type = rep(NA, 2*N))
  plot_df$m6A_ratio[1:N] <- assays(se)$m6A[indx,]/assays(se)$Total[indx,]
  plot_df$type[1:N] <- "estimated from data"
  p <- sum(m6A[indx])/sum(Total[indx])
  plot_df$m6A_ratio[(N+1):(2*N)] <- rbinom(n = N, size = Total[indx], prob = p)/Total[indx]
  plot_df$type[(N+1):(2*N)] <- "simulated by fitted model"
  qs <- 1-rev(quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], probs = seq(0,1,length.out = 250)))
  qtl_emp <- quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], qs)
  qtl_fit <- quantile(plot_df$m6A_ratio[plot_df$type != "estimated from data"], qs)
  D <- 1
  BIC <- D*log(N) - 2*sum(dbinom(m6A[indx], Total[indx], p, log = TRUE))
  return( list(plot_df_qtl = data.frame(qtl_emp = qtl_emp,
                                        qtl_fit = qtl_fit),
               plot_df_hist = plot_df,
               BIC = BIC) )
}

plot_hist <- function(plot_df, model_name){
  ggplot(plot_df, aes(x=m6A_ratio,fill=type)) +
    geom_histogram(colour = "dark blue", binwidth = 0.025,size = 0.2) +
    facet_wrap(vars(type)) + theme_bw() + scale_fill_brewer(palette = "Accent") + labs(x = "m6A Ratio", y = "Count")
  ggsave(paste0(model_name, "_hist.pdf"), width = 7.5, height = 2.5)
}

extract_bbFit <- function(m6A, Total, total_threshold = 30){
  #require(SummarizedExperiment)
  se <- SummarizedExperiment(assays = list(m6A = cbind(m6A), Total = cbind(Total)), checkDimnames = FALSE)
  indx <- as.vector(assays(se)$Total >= total_threshold)
  N <- sum(indx)
  plot_df <- data.frame(m6A_ratio = rep(NA, 2*N), type = rep(NA, 2*N))
  plot_df$m6A_ratio[1:N] <- assays(se)$m6A[indx,]/assays(se)$Total[indx,]
  plot_df$type[1:N] <- "estimated from data"
  bbfit <- beta_mom(m6A[indx], Total[indx], w = rep(1, sum(indx)))
  bbfit <- BB_WMLE(m6A[indx], Total[indx], w = rep(1, sum(indx)), a = bbfit[1], b = bbfit[2])
  plot_df$m6A_ratio[(N+1):(2*N)] <- rbbinom(n = N, size = Total[indx], alpha = bbfit[1], beta = bbfit[2])/Total[indx]
  plot_df$type[(N+1):(2*N)] <- "simulated by fitted model"
  qs <- 1-rev(quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], probs = seq(0,1,length.out = 250)))
  qtl_emp <- quantile(plot_df$m6A_ratio[plot_df$type == "estimated from data"], qs)
  qtl_fit <- quantile(plot_df$m6A_ratio[plot_df$type != "estimated from data"], qs)
  D <- 2
  BIC <- D*log(N) - 2*sum(log(dbbinom(m6A[indx], Total[indx], bbfit[1], bbfit[2])))
  return( list(plot_df_qtl = data.frame(qtl_emp = qtl_emp,
                                        qtl_fit = qtl_fit),
               plot_df_hist = plot_df,
               BIC = BIC) )
}

#' Compare Goodness of Fits for m6A Modification Models
#'
#' This function compares the goodness of fit between several model types for m6A modification data. It generates and saves two types of plots: histograms of beta/m6A level values and quantile-quantile (QQ) plots to assess model fits.
#'
#' @param m6A A numeric vector of counts for m6A modifications.
#' @param Total A numeric vector of total counts.
#' @param cov_threshold A numeric threshold to filter sites based on coverage, only sites with total read count >= `cov_threshold` will be used in evaluation and model fitting. Default is 30.
#' @param bbmix_size An integer specifying the maximum number of sites for the bbmix model. This parameter is for efficiency considerations and only affects the "bbmix" method. Default is NULL.
#'
#' @return The function outputs two figures in the current working directory:
#' 1. A histogram plot ("compareGof_hist.pdf") illustrating the degree of similarity between empirical data and simulations from fitted models.
#' 2. A QQ plot ("compareGof_qq.pdf") comparing the goodness of fits of different models.
#' 3. A CSV table ("compareGof_BIC.csv") comparing the Bayesian Information Criterion (BIC) of different models.
#'
#' @details
#' The function generates histograms for each model type, including binomial, binomial mixture, binomial uniform mixture, and beta-binomial mixture.
#' The histograms compare the observed m6A ratios to the model predictions.
#' The QQ plot illustrates the fit of the quantiles from the empirical data against the fitted models.
#' The CSV table contains the BIC values of the five fitted models.
#'
#' @import ggplot2
#' @export
#' @examples
#' library(SummarizedExperiment)
#' m6A_se <- readRDS(system.file("extdata", "example_se.rds", package="OmixM6A"))
#' m6A_counts <- assays(m6A_se)$m6A[,1]
#' total_counts <- assays(m6A_se)$Total[,1]
#' set.seed(123)
#' compareGoodnessOfFits(m6A_counts, total_counts)
compareGoodnessOfFits <- function(m6A, Total, cov_threshold = 30, bbmix_size = NULL){
  #require(ggplot2)
  binomFit_ls <- extract_binomialFit(m6A, Total,  cov_threshold)
  bbFit_ls <- extract_bbFit(m6A, Total,  cov_threshold)
  bmixFit_ls <- extract_bmixFit(m6A, Total,  cov_threshold)
  bumixFit_ls <- extract_bumixFit(m6A, Total,  cov_threshold)
  bbmixFit_ls <- extract_bbmixFit(m6A, Total,  cov_threshold)
  binomFit_ls$plot_df_hist$type[binomFit_ls$plot_df_hist$type != "estimated from data"] <- "binomial"
  bbFit_ls$plot_df_hist$type[bbFit_ls$plot_df_hist$type != "estimated from data"] <- "beta-binomial"
  bmixFit_ls$plot_df_hist$type[bmixFit_ls$plot_df_hist$type != "estimated from data"] <- "binomial mixture"
  bumixFit_ls$plot_df_hist$type[bumixFit_ls$plot_df_hist$type != "estimated from data"] <- "binomial uniform mixture"
  bbmixFit_ls$plot_df_hist$type[bbmixFit_ls$plot_df_hist$type != "estimated from data"] <- "beta-binomial mixture"
  plot_df_hist <- rbind(binomFit_ls$plot_df_hist,
                        bbFit_ls$plot_df_hist,
                        bmixFit_ls$plot_df_hist,
                        bumixFit_ls$plot_df_hist)
  plot_df_hist <- plot_df_hist[plot_df_hist$type != "estimated from data",]
  plot_df_hist <- rbind(plot_df_hist, bbmixFit_ls$plot_df_hist)
  plot_df_hist$type <- factor(plot_df_hist$type, levels = c("estimated from data", "beta-binomial mixture", "beta-binomial", "binomial", "binomial mixture", "binomial uniform mixture"))
  ggplot(plot_df_hist, aes(x=m6A_ratio,fill=type)) +
    geom_histogram(colour = "dark blue", binwidth = 0.025,size = 0.2) +
    facet_wrap(vars(type)) + theme_bw() + scale_fill_brewer(palette = "Accent") + labs(x = "m6A Ratio", y = "Count")
  ggsave("compareGof_hist.pdf", width = 7.8, height = 3.5)

  binomFit_ls$plot_df_qtl$model <- "binomial"
  bmixFit_ls$plot_df_qtl$model <- "binomial mixture"
  bumixFit_ls$plot_df_qtl$model <- "binomial uniform mixture"
  bbFit_ls$plot_df_qtl$model <- "beta-binomial"
  bbmixFit_ls$plot_df_qtl$model <- "beta-binomial mixture"

  plot_df_qq <- rbind(binomFit_ls$plot_df_qtl,
                      bbFit_ls$plot_df_qtl,
                      bmixFit_ls$plot_df_qtl,
                      bumixFit_ls$plot_df_qtl,
                      bbmixFit_ls$plot_df_qtl)
  plot_df_qq$model <- factor(plot_df_qq$model, levels = c("binomial", "binomial mixture", "binomial uniform mixture", "beta-binomial", "beta-binomial mixture"))
  ggplot(plot_df_qq, aes(x=qtl_fit, y=qtl_emp, colour = model, shape = model)) +
    geom_point(size = 0.8) +
    theme_classic() +
    labs(x = "Fitted Quantiles", y = "Emperical Quantiles") + xlim(0,1) + ylim(0,1) + scale_colour_brewer(palette = "Dark2")
  ggsave("compareGof_qq.pdf", width = 5, height = 2.5)

  BIC_tbl <- data.frame(model = c("binomial", "binomial mixture", "binomial uniform mixture", "beta-binomial", "beta-binomial mixture"),
                        BIC = c(binomFit_ls$BIC, bmixFit_ls$BIC, bumixFit_ls$BIC, bbFit_ls$BIC, bbmixFit_ls$BIC))
  write.csv(BIC_tbl, "compareGof_BIC.csv", row.names = FALSE)
}
