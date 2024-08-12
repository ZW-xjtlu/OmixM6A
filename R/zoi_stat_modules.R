# The function fits zero & one inflated 2 component binomial mixture and return the fitted parameters together with posterior for foreground
fit_zoibmix <- function(m6A_vec, Total_vec, cov_threshold){
  indx_zero <- Total_vec == 0
  indx_sub <- Total_vec[!indx_zero] >= cov_threshold
  m6A <- m6A_vec[!indx_zero][indx_sub]
  Ns <- Total_vec[!indx_zero][indx_sub]

  #Fit a frequentist EM algorithm to estimate parameters from data (m6A and Ns)
  # Use closer initial guesses to the true values
  bg_prop <- 0.85
  fg_prop <- 0.10
  zero_prop <- 0.04
  one_prop <- 0.01
  p_bg <- 0.04
  p_fg <- 0.45

  # Initialize
  tolerance <- 1e-5
  max_iter <- 100

  # Track previous parameter values
  prev_p_bg <- 0
  prev_p_fg <- 0
  prev_bg_prop <- 0
  prev_fg_prop <- 0
  prev_zero_prop <- 0
  prev_one_prop <- 0
  m <- length(m6A)

  for(i in 1:max_iter){
    # E-step in log-space
    log_prob_zero <- ifelse(m6A == 0, log(zero_prop), -Inf)
    log_prob_one <- ifelse(m6A == Ns, log(one_prop), -Inf)
    log_prob_fg <- dbinom(m6A, size = Ns, prob = p_fg, log = TRUE) + log(fg_prop)
    log_prob_bg <- dbinom(m6A, size = Ns, prob = p_bg, log = TRUE) + log(bg_prop)
    max_log_prob <- pmax(log_prob_fg, log_prob_bg, log_prob_zero, log_prob_one)

    # Calculate responsibilities while avoiding underflow
    denum <- exp(log_prob_fg - max_log_prob) + exp(log_prob_bg - max_log_prob) +  exp(log_prob_zero - max_log_prob) + exp(log_prob_one - max_log_prob)
    resp_fg <- exp(log_prob_fg - max_log_prob) / denum
    resp_bg <- exp(log_prob_bg - max_log_prob) / denum
    resp_zero <- exp(log_prob_zero - max_log_prob) / denum
    resp_one <- exp(log_prob_one - max_log_prob) / denum

    # M-step weighted MLE of binomial parameters
    p_bg <- sum((resp_bg) * m6A) / sum((resp_bg) * Ns)
    p_fg <- sum(resp_fg * m6A) / sum(resp_fg * Ns)

    # Update the proportions based on responsibilities
    zero_prop <- sum(resp_zero) / m
    one_prop <- sum(resp_one) / m
    bg_prop <- sum(resp_bg) / m
    fg_prop <- sum(resp_fg) / m

    # Check for convergence
    if (abs(p_bg - prev_p_bg) < tolerance &&
        abs(p_fg - prev_p_fg) < tolerance &&
        abs(fg_prop - prev_fg_prop) < tolerance &&
        abs(bg_prop - prev_bg_prop) < tolerance &&
        abs(zero_prop - prev_zero_prop) < tolerance &&
        abs(one_prop - prev_one_prop) < tolerance) {
      cat("Converged in", i, "iterations.\n")
      break
    }

    # Update previous values
    prev_p_bg <- p_bg
    prev_p_fg <- p_fg
    prev_bg_prop <- bg_prop
    prev_fg_prop <- fg_prop
    prev_zero_prop <- zero_prop
    prev_one_prop <- one_prop
  }

  #Check and correct switching between foreground & background components
  if (p_bg > p_fg) {
    tmp <- p_bg
    p_bg <- p_fg
    p_fg <- tmp
    tmp <- resp_bg
    resp_bg <- resp_fg
    resp_fg <- tmp
  }

  #The returned posterior is the probability of foreground component + one component
  resp_return <- rep(NA, length(Total_vec))
  resp_return[!indx_zero] <- resp_fg + resp_one
  prob_bn_fg <- rep(NA, length(Total_vec))
  prob_bn_fg[!indx_zero] <- resp_fg
  prob_bn_bg <- rep(NA, length(Total_vec))
  prob_bn_bg[!indx_zero] <- resp_bg

  #Caculaate p-value of zero-inflated binomial
  pvalue <- rep(NA, length(Total_vec))
  pi <- resp_bg/(resp_bg+resp_zero)
  pi[is.na(pi)] <- 0
  pvalue[!indx_zero] <- pbinom(m6A_vec[!indx_zero]-1, Total_vec[!indx_zero],
                               p_bg,
                               lower.tail = FALSE) * pi
  pvalue[!indx_zero][m6A_vec[!indx_zero] == 0] <- 1
  beta <- rep(NA, length(Total_vec))
  beta[!indx_zero] <- m6A_vec[!indx_zero]/Total_vec[!indx_zero]
  output_lst <- list(prob_fg = resp_return,
                     prob_bn_fg = prob_bn_fg,
                     prob_bn_bg = prob_bn_bg,
                     pvalue = pvalue,
                     beta = beta,
                     para = list(zero_proportion = zero_prop,
                                 bg_proportion = bg_prop,
                                 fg_proportion = fg_prop,
                                 one_proportion = one_prop,
                                 p_m6A_bg = p_bg,
                                 p_m6A_fg = p_fg))
  return(output_lst)
}

# The function fits zero & one inflated 2 components beta binomial mixture and return the fitted parameters together with posterior for foreground
fit_zoibbmix <- function(m6A_vec, Total_vec, subsize = NULL, cov_threshold){
  #require(extraDistr)
  indx_zero <- Total_vec == 0
  #Subset data
  if(!is.null(subsize)) {
    indx_sub <- Total_vec[!indx_zero] >= cov_threshold
    indx_sub2 <- sample.int(sum(indx_sub), min(sum(indx_sub), subsize))
    m6A <- m6A_vec[!indx_zero][indx_sub][indx_sub2]
    Ns <- Total_vec[!indx_zero][indx_sub][indx_sub2]
  }else{
    indx_sub <- Total_vec[!indx_zero] >= cov_threshold
    m6A <- m6A_vec[!indx_zero][indx_sub]
    Ns <- Total_vec[!indx_zero][indx_sub]
  }
  #Stabilize huge count if any
  indx_huge <- Ns > 100000
  if(any(indx_huge)){
    ratio_huge <- m6A[indx_huge]/Ns[indx_huge]
    Ns[indx_huge] <- 100000
    m6A[indx_huge] <- round(ratio_huge*Ns[indx_huge])
  }

  # Obtain initial estimates using b_mix
  zibmix_fit <- fit_zoibmix( m6A, Ns, cov_threshold )
  #Fit an EM algorithm to estimate parameters from data (m6A and Ns)
  # Use closer initial guesses to the true values
  zero_prop <- zibmix_fit$para$zero_proportion
  bg_prop <- zibmix_fit$para$bg_proportion
  fg_prop <- zibmix_fit$para$fg_proportion
  one_prop <- zibmix_fit$para$one_proportion
  # Set initial values with MoM
  mom_fit_fg <- beta_mom(m6A, Ns, zibmix_fit$prob_fg, 30)
  wmle_fg <- BB_WMLE(m6A, Ns, zibmix_fit$prob_bn_fg, mom_fit_fg[1], mom_fit_fg[2])
  mom_fit_bg <- beta_mom(m6A, Ns, 1-zibmix_fit$prob_fg, 30)
  wmle_bg <- BB_WMLE(m6A, Ns, zibmix_fit$prob_bn_bg, mom_fit_bg[1], mom_fit_bg[2])

  # Detect broken MLE results
  if(any(c(wmle_bg, wmle_fg) < 0) | any(c(wmle_bg, wmle_fg) > 1000)) {
    warning("MLE has failed convergence, MoM is used instead in one EM iteration.")
    wmle_fg <- mom_fit_fg
    wmle_bg <- mom_fit_bg
  }

  alpha_fg <- wmle_fg[1]
  beta_fg <- wmle_fg[2]
  alpha_bg <- wmle_bg[1]
  beta_bg <- wmle_bg[2]

  # Initialize
  tolerance <- 1e-3
  max_iter <- 1000

  # Track previous parameter values
  prev_alpha_fg <- 0
  prev_beta_fg <- 0
  prev_alpha_bg <- 0
  prev_beta_bg <- 0
  prev_zero_prop <- 0
  prev_bg_prop <- 0
  prev_fg_prop <- 0
  prev_one_prop <- 0
  m <- length(m6A)

  for(i in 1:max_iter){
    cat("Run EM iteration", i, "...\n")
    cat("Current alpha fg: ", alpha_fg, "...\n")
    cat("Current beta fg: ", beta_fg, "...\n")
    cat("Current alpha bg: ", alpha_bg, "...\n")
    cat("Current beta bg: ", beta_bg, "...\n")
    cat("Current zero proportion: ", zero_prop, "...\n")
    cat("Current fg proportion: ", fg_prop, "...\n")
    cat("Current bg proportion: ", bg_prop, "...\n")
    cat("Current one proportion: ", one_prop, "...\n")

    # E-step in log-space
    #log_prob_fg <- dbinom(m6A, size = Ns, prob = p_fg, log = TRUE) + log(fg_prop)
    #log_prob_bg <- dbinom(m6A, size = Ns, prob = p_bg, log = TRUE) + log(bg_prop)
    log_prob_zero <- ifelse(m6A == 0, log(zero_prop), -Inf)
    log_prob_one <- ifelse(m6A == Ns, log(one_prop), -Inf)
    log_prob_fg <- dbbinom(m6A, Ns, alpha_fg, beta_fg, log = TRUE) + log(fg_prop)
    log_prob_bg <- dbbinom(m6A, Ns, alpha_bg, beta_bg, log = TRUE) + log(bg_prop)
    max_log_prob <- pmax(log_prob_fg, log_prob_bg, log_prob_zero, log_prob_one)


    # Calculate responsibilities while avoiding underflow
    denum <- exp(log_prob_fg - max_log_prob) + exp(log_prob_bg - max_log_prob) +  exp(log_prob_zero - max_log_prob) + exp(log_prob_one - max_log_prob)
    resp_fg <- exp(log_prob_fg - max_log_prob) / denum
    resp_bg <- exp(log_prob_bg - max_log_prob) / denum
    resp_zero <- exp(log_prob_zero - max_log_prob) / denum
    resp_one <- exp(log_prob_one - max_log_prob) / denum

    # M-step
    wmle_fg <- BB_WMLE(m6A, Ns, resp_fg, alpha_fg, beta_fg)
    wmle_bg <- BB_WMLE(m6A, Ns, resp_bg, alpha_bg, beta_bg)

    # Detect broken MLE results
    if(any(c(wmle_bg, wmle_fg) < 0) | any(c(wmle_bg, wmle_fg) > 1000)) {
      warning("MLE has failed convergence, MoM is used instead in one EM iteration.")
      wmle_fg <- beta_mom(m6A, Ns, resp_fg, 30)
      wmle_bg <- beta_mom(m6A, Ns, resp_bg, 30)
    }

    alpha_fg <- wmle_fg[1]
    beta_fg <- wmle_fg[2]
    alpha_bg <- wmle_bg[1]
    beta_bg <- wmle_bg[2]

    # Update the proportions based on responsibilities
    zero_prop <- sum(resp_zero) / m
    one_prop <- sum(resp_one) / m
    bg_prop <- sum(resp_bg) / m
    fg_prop <- sum(resp_fg) / m

    # Check for convergence
    if (abs(alpha_bg - prev_alpha_bg) < tolerance &&
        abs(beta_bg - prev_beta_bg) < tolerance &&
        abs(alpha_fg - prev_alpha_fg) < tolerance &&
        abs(beta_fg - prev_beta_fg) < tolerance &&
        abs(fg_prop - prev_fg_prop) < tolerance &&
        abs(bg_prop - prev_bg_prop) < tolerance &&
        abs(zero_prop - prev_zero_prop) < tolerance &&
        abs(one_prop - prev_one_prop) < tolerance) {
      cat("Converged in", i, "iterations.\n")
      break
    }

    # Update previous values
    prev_alpha_bg <- alpha_bg
    prev_alpha_fg <- alpha_fg
    prev_beta_bg <- beta_bg
    prev_beta_fg <- beta_fg
    prev_bg_prop <- bg_prop
    prev_fg_prop <- fg_prop
    prev_zero_prop <- zero_prop
    prev_one_prop <- one_prop
  }

  #Check consistency of foreground & background components
  if ((alpha_bg/(alpha_bg+beta_bg)) > (alpha_fg/(alpha_fg+beta_fg))) {
    tmp <- alpha_bg
    alpha_bg <- alpha_fg
    alpha_fg <- tmp
    tmp <- beta_bg
    beta_bg <- beta_fg
    beta_fg <- tmp
  }

  #Update responsibilities on all non-zero sites
  m6A <- m6A_vec[!indx_zero]
  Ns <- Total_vec[!indx_zero]
  log_prob_zero <- ifelse(m6A == 0, log(zero_prop), -Inf)
  log_prob_one <- ifelse(m6A == Ns, log(one_prop), -Inf)
  log_prob_fg <- dbbinom(m6A, Ns, alpha_fg, beta_fg, log = TRUE) + log(fg_prop)
  log_prob_bg <- dbbinom(m6A, Ns, alpha_bg, beta_bg, log = TRUE) + log(bg_prop)
  max_log_prob <- pmax(log_prob_fg, log_prob_bg, log_prob_zero, log_prob_one)
  denum <- exp(log_prob_fg - max_log_prob) + exp(log_prob_bg - max_log_prob) +  exp(log_prob_zero - max_log_prob) + exp(log_prob_one - max_log_prob)
  resp_fg <- exp(log_prob_fg - max_log_prob) / denum
  resp_bg <- exp(log_prob_bg - max_log_prob) / denum
  resp_zero <- exp(log_prob_zero - max_log_prob) / denum
  resp_one <- exp(log_prob_one - max_log_prob) / denum
  resp_return <- rep(NA, length(Total_vec))
  resp_return[!indx_zero] <- resp_fg + resp_one

  #Caculaate p-value of zero-inflated beta-binomial
  pvalue <- rep(NA, length(Total_vec))
  pi <- resp_bg/(resp_bg+resp_zero)
  pi[is.na(pi)] <- 0
  pvalue[!indx_zero] <- pbbinom(m6A_vec[!indx_zero]-1, Total_vec[!indx_zero],
                               alpha_bg, beta_bg,
                               lower.tail = FALSE) * pi
  pvalue[!indx_zero][m6A_vec[!indx_zero] == 0] <- 1
  beta <- rep(NA, length(Total_vec))
  beta[!indx_zero] <- m6A_vec[!indx_zero]/Total_vec[!indx_zero]
  output_lst <- list(prob_fg = resp_return,
                     pvalue = pvalue,
                     beta = beta,
                     para = list(zero_proportion = zero_prop,
                                 bg_proportion = bg_prop,
                                 fg_proportion = fg_prop,
                                 one_proportion = one_prop,
                                 alpha_m6A_bg = alpha_bg,
                                 beta_m6A_bg = beta_bg,
                                 alpha_m6A_fg = alpha_fg,
                                 beta_m6A_fg = beta_fg))
  return(output_lst)
}
