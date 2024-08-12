# Weighted method of moments estimates of beta-binomial parameters
beta_mom <- function(x, n, w, min_data_points = 10, desired_quantile = 0.75){
  # Determine the threshold for filtering based on quantiles if there are too few data points above 30
  n_filter <- if(sum(n > 30) < min_data_points) {
    quantile(n, probs = desired_quantile, na.rm = TRUE)
  } else {
    30
  }
  indx <- n >= n_filter
  sample_prop <- x[indx]/n[indx]
  # Calculate the weighted mean
  weighted_mean <- sum(w[indx] * sample_prop) / sum(w[indx])

  # Calculate the weighted variance
  weighted_var <- sum(w[indx] * (sample_prop - weighted_mean)^2) / sum(w[indx])

  # Method of Moments estimates for alpha and beta
  alpha_mom <- weighted_mean * ((weighted_mean * (1 - weighted_mean) / weighted_var) - 1)
  beta_mom <- (1 - weighted_mean) * ((weighted_mean * (1 - weighted_mean) / weighted_var) - 1)

  return(c(alpha_mom, beta_mom))
}

BB_WMLE <- function(x, n, w, a, b){
  # Log-likelihood partial derivatives with respect to alpha, beta
  ll.a <- function(a, b, w) sum(w * (digamma(a + b) - digamma(a) + digamma(a + x) - digamma(a + b + n)))
  ll.b <- function(a, b, w) sum(w * (digamma(a + b) - digamma(b) + digamma(b + n - x) - digamma(a + b + n)))
  # Second derivatives of the log-likelihood
  ll.aa <- function(a, b, w) sum(w * (trigamma(a + b) - trigamma(a) + trigamma(a + x) - trigamma(a + b + n)))
  ll.bb <- function(a, b, w) sum(w * (trigamma(a + b) - trigamma(b) + trigamma(b + n - x) - trigamma(a + b + n)))
  ll.ab <- function(a, b, w) sum(w * (trigamma(a + b) - trigamma(a + b + n)))
  # Optimization loop
  a2 <- a; b2 <- b # Assuming 'a' and 'b' are initial estimates for the parameters alpha and beta
  for(i in 1:50) {
    parm <- matrix(c(a2, b2), nrow = 2)
    f <- matrix(c(ll.a(a2, b2, w), ll.b(a2, b2, w)), nrow = 2)
    J <- ginv(matrix(c(ll.aa(a2, b2, w), ll.ab(a2, b2, w),
                        ll.ab(a2, b2, w), ll.bb(a2, b2, w)),
                      nrow = 2, ncol = 2))
    f2 <- parm - J %*% f
    a2 <- f2[1, ]; b2 <- f2[2, ]
    if(max(abs(c(parm[1, ] - a2, parm[2, ] - b2))) < 1e-6) break
  }
  return(c(a2, b2))
}

# The function fits 2 components beta binomial mixture and return the fitted parameters together with posterior for foreground
fit_bbmix <- function(m6A_vec, Total_vec, subsize = NULL, cov_threshold){
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
  bmix_fit <- fit_bmix( m6A, Ns, cov_threshold)
  #Fit an EM algorithm to estimate parameters from data (m6A and Ns)
  # Use closer initial guesses to the true values
  bg_prop <- bmix_fit$para$bg_proportion
  fg_prop <- bmix_fit$para$fg_proportion
  # Set initial values with MoM
  mom_fit_fg <- beta_mom(m6A, Ns, bmix_fit$prob_fg, 30)
  wmle_fg <- BB_WMLE(m6A, Ns, bmix_fit$prob_fg, mom_fit_fg[1], mom_fit_fg[2])
  mom_fit_bg <- beta_mom(m6A, Ns, 1-bmix_fit$prob_fg, 30)
  wmle_bg <- BB_WMLE(m6A, Ns, 1-bmix_fit$prob_fg, mom_fit_bg[1], mom_fit_bg[2])

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
  prev_bg_prop <- 0
  prev_fg_prop <- 0
  m <- length(m6A)

  for(i in 1:max_iter){
    cat("Run EM iteration", i, "...\n")
    cat("Current alpha fg: ", alpha_fg, "...\n")
    cat("Current beta fg: ", beta_fg, "...\n")
    cat("Current alpha bg: ", alpha_bg, "...\n")
    cat("Current beta bg: ", beta_bg, "...\n")
    cat("Current fg proportion: ", fg_prop, "...\n")
    cat("Current bg proportion: ", bg_prop, "...\n")

    # E-step in log-space
    #log_prob_fg <- dbinom(m6A, size = Ns, prob = p_fg, log = TRUE) + log(fg_prop)
    #log_prob_bg <- dbinom(m6A, size = Ns, prob = p_bg, log = TRUE) + log(bg_prop)
    log_prob_fg <- dbbinom(m6A, Ns, alpha_fg, beta_fg, log = TRUE) + log(fg_prop)
    log_prob_bg <- dbbinom(m6A, Ns, alpha_bg, beta_bg, log = TRUE) + log(bg_prop)
    max_log_prob <- pmax(log_prob_fg, log_prob_bg)


    # Calculate responsibilities while avoiding underflow
    responsibilities <- exp(log_prob_fg - max_log_prob) / (exp(log_prob_fg - max_log_prob) + exp(log_prob_bg - max_log_prob))

    # M-step
      wmle_fg <- BB_WMLE(m6A, Ns, responsibilities, alpha_fg, beta_fg)
      wmle_bg <- BB_WMLE(m6A, Ns, 1-responsibilities, alpha_bg, beta_bg)

    # Detect broken MLE results
    if(any(c(wmle_bg, wmle_fg) < 0) | any(c(wmle_bg, wmle_fg) > 1000)) {
      warning("MLE has failed convergence, MoM is used instead in one EM iteration.")
      wmle_fg <- beta_mom(m6A, Ns, responsibilities, 30)
      wmle_bg <- beta_mom(m6A, Ns, 1-responsibilities, 30)
    }

    alpha_fg <- wmle_fg[1]
    beta_fg <- wmle_fg[2]
    alpha_bg <- wmle_bg[1]
    beta_bg <- wmle_bg[2]

    # Update the proportions based on responsibilities
    bg_prop <- sum(1 - responsibilities) / m
    fg_prop <- sum(responsibilities) / m

    # Check for convergence
    if (abs(alpha_bg - prev_alpha_bg) < tolerance &&
        abs(beta_bg - prev_beta_bg) < tolerance &&
        abs(alpha_fg - prev_alpha_fg) < tolerance &&
        abs(beta_fg - prev_beta_fg) < tolerance &&
        abs(fg_prop - prev_fg_prop) < tolerance &&
        abs(bg_prop - prev_bg_prop) < tolerance) {
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
  log_prob_fg <- dbbinom(m6A_vec[!indx_zero], Total_vec[!indx_zero], alpha_fg, beta_fg, log = TRUE) + log(fg_prop)
  log_prob_bg <- dbbinom(m6A_vec[!indx_zero], Total_vec[!indx_zero], alpha_bg, beta_bg, log = TRUE) + log(bg_prop)
  max_log_prob <- pmax(log_prob_fg, log_prob_bg)
  responsibilities <- exp(log_prob_fg - max_log_prob) / (exp(log_prob_fg - max_log_prob) + exp(log_prob_bg - max_log_prob))
  resp_return <- rep(NA, length(Total_vec))
  resp_return[!indx_zero] <- responsibilities

  pvalue <- rep(NA, length(Total_vec))
  pvalue[!indx_zero] <- pbbinom(m6A_vec[!indx_zero]-1, Total_vec[!indx_zero],
                                alpha_bg, beta_bg,
                                lower.tail = FALSE)

  beta <- rep(NA, length(Total_vec))
  beta[!indx_zero] <- m6A_vec[!indx_zero]/Total_vec[!indx_zero]
  output_lst <- list(prob_fg = resp_return,
                     pvalue = pvalue,
                     beta = beta,
                     para = list(bg_proportion = bg_prop,
                                 fg_proportion = fg_prop,
                                 alpha_m6A_bg = alpha_bg,
                                 beta_m6A_bg = beta_bg,
                                 alpha_m6A_fg = alpha_fg,
                                 beta_m6A_fg = beta_fg))
  return(output_lst)
}

# The function fits 2 component binomial mixture and return the fitted parameters together with posterior for foreground
fit_bmix <- function(m6A_vec, Total_vec, cov_threshold){
  indx_zero <- Total_vec == 0
  indx_sub <- Total_vec[!indx_zero] >= cov_threshold
  m6A <- m6A_vec[!indx_zero][indx_sub]
  Ns <- Total_vec[!indx_zero][indx_sub]

  #Fit a frequentist EM algorithm to estimate parameters from data (m6A and Ns)
  # Use closer initial guesses to the true values
  bg_prop <- 0.85
  fg_prop <- 0.15
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
  m <- length(m6A)

  for(i in 1:max_iter){
    # E-step in log-space
    log_prob_fg <- dbinom(m6A, size = Ns, prob = p_fg, log = TRUE) + log(fg_prop)
    log_prob_bg <- dbinom(m6A, size = Ns, prob = p_bg, log = TRUE) + log(bg_prop)
    max_log_prob <- pmax(log_prob_fg, log_prob_bg)

    # Calculate responsibilities while avoiding underflow
    responsibilities <- exp(log_prob_fg - max_log_prob) / (exp(log_prob_fg - max_log_prob) + exp(log_prob_bg - max_log_prob))

    # M-step
    p_bg <- sum((1 - responsibilities) * m6A) / sum((1 - responsibilities) * Ns)
    p_fg <- sum(responsibilities * m6A) / sum(responsibilities * Ns)

    # Update the proportions based on responsibilities
    bg_prop <- sum(1 - responsibilities) / m
    fg_prop <- sum(responsibilities) / m

    # Check for convergence
    if (abs(p_bg - prev_p_bg) < tolerance &&
        abs(p_fg - prev_p_fg) < tolerance &&
        abs(fg_prop - prev_fg_prop) < tolerance &&
        abs(bg_prop - prev_bg_prop) < tolerance) {
      cat("Converged in", i, "iterations.\n")
      break
    }

    # Update previous values
    prev_p_bg <- p_bg
    prev_p_fg <- p_fg
    prev_bg_prop <- bg_prop
    prev_fg_prop <- fg_prop
  }

  resp_return <- rep(NA, length(Total_vec))
  resp_return[!indx_zero] <- responsibilities
  pvalue <- rep(NA, length(Total_vec))
  pvalue[!indx_zero] <- pbinom(m6A_vec[!indx_zero]-1, Total_vec[!indx_zero],
                               p_bg,
                               lower.tail = FALSE)
  beta <- rep(NA, length(Total_vec))
  beta[!indx_zero] <- m6A_vec[!indx_zero]/Total_vec[!indx_zero]
  output_lst <- list(prob_fg = resp_return,
                     pvalue = pvalue,
                     beta = beta,
                     para = list(bg_proportion = bg_prop,
                                 fg_proportion = fg_prop,
                                 p_m6A_bg = p_bg,
                                 p_m6A_fg = p_fg))
  return(output_lst)
}


# The function fits 2 component mixture model of binomial and beta(1, 1) / uniform binomial distribution
fit_bumix <- function(m6A_vec, Total_vec, cov_threshold){
  #require(extraDistr)
  indx_zero <- Total_vec == 0
  indx_sub <- Total_vec[!indx_zero] >= cov_threshold
  m6A <- m6A_vec[!indx_zero][indx_sub]
  Ns <- Total_vec[!indx_zero][indx_sub]

  #Fit a frequentist EM algorithm to estimate parameters from data (m6A and Ns)
  # Use closer initial guesses to the true values
  bg_prop <- 0.85
  fg_prop <- 0.15
  p_bg <- 0.04

  # Initialize
  tolerance <- 1e-5
  max_iter <- 100

  # Track previous parameter values
  prev_p_bg <- 0
  prev_bg_prop <- 0
  prev_fg_prop <- 0
  m <- length(m6A)

  for(i in 1:max_iter){
    # E-step in log-space
    log_prob_fg <- dbbinom(m6A, Ns, 1, 1, log = TRUE) + log(fg_prop)
    log_prob_bg <- dbinom(m6A, size = Ns, prob = p_bg, log = TRUE) + log(bg_prop)
    max_log_prob <- pmax(log_prob_fg, log_prob_bg)

    # Calculate responsibilities while avoiding underflow
    responsibilities <- exp(log_prob_fg - max_log_prob) / (exp(log_prob_fg - max_log_prob) + exp(log_prob_bg - max_log_prob))

    # M-step
    p_bg <- sum((1 - responsibilities) * m6A) / sum((1 - responsibilities) * Ns)

    # Update the proportions based on responsibilities
    bg_prop <- sum(1 - responsibilities) / m
    fg_prop <- sum(responsibilities) / m

    # Check for convergence
    if ((abs(p_bg - prev_p_bg) < tolerance &&
        abs(fg_prop - prev_fg_prop) < tolerance &&
        abs(bg_prop - prev_bg_prop) < tolerance)|(fg_prop == 1)) {
      cat("Converged in", i, "iterations.\n")
      break
    }

    # Update previous values
    prev_p_bg <- p_bg
    prev_bg_prop <- bg_prop
    prev_fg_prop <- fg_prop
  }

  resp_return <- rep(NA, length(Total_vec))
  resp_return[!indx_zero] <- responsibilities
  pvalue <- rep(NA, length(Total_vec))
  pvalue[!indx_zero] <- pbinom(m6A_vec[!indx_zero]-1, Total_vec[!indx_zero],
                               p_bg,
                               lower.tail = FALSE)
  beta <- rep(NA, length(Total_vec))
  beta[!indx_zero] <- m6A_vec[!indx_zero]/Total_vec[!indx_zero]
  output_lst <- list(prob_fg = resp_return,
                     pvalue = pvalue,
                     beta = beta,
                     para = list(bg_proportion = bg_prop,
                                 fg_proportion = fg_prop,
                                 p_m6A_bg = p_bg))
  return(output_lst)
}

