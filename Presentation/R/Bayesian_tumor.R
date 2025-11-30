#' Bayesian Inference for Brain Tumor Growth Parameters with Parallel Chains
#'
#' This function provides Bayesian parameter estimation for reaction-diffusion
#' models of brain tumor growth using adaptive MCMC sampling with parallel chains.
#'
#' @param noise_std Observation noise standard deviation. 
#' @param t_obs Observation time in days 
#' @param n_chains Number of MCMC chains 
#' @param n_iter MCMC iterations per chain 
#' @param burnin Burn-in period per chain 
#' @param grid_size Number of spatial grid points
#' @param domain_length Domain length in mm 
#' @param dt Time step for PDE solver 
#' @param true_D True diffusion coefficient 
#' @param true_rho True proliferation rate 
#' @param use_parallel Use parallel computation for multiple chains 
#' @param seed Random seed
#'
#' @return A list containing MCMC results, diagnostics, and performance metrics
#'
#' @export
bayesian_tumor_growth <- function(noise_std = 0.05, t_obs = 30, n_chains = 4, n_iter = 5000,
                                  burnin = 1000, grid_size = 100, domain_length = 50,
                                  dt = 0.01, true_D = 1.0, true_rho = 0.18, use_parallel = TRUE,
                                  seed = NULL) {
  # Check input
  if (!is.numeric(noise_std) || noise_std <= 0) {
    stop("noise_std must be a positive numeric value.")
  }
  
  if (!is.numeric(t_obs) || t_obs <= 0) {
    stop("t_obs must be a positive numeric value.")
  }
  
  if (!is.numeric(n_chains) || n_chains <= 0 || n_chains != as.integer(n_chains)) {
    stop("n_chains must be a positive integer.")
  }
  
  if (!is.numeric(n_iter) || n_iter <= 0 || n_iter != as.integer(n_iter)) {
    stop("n_iter must be a positive integer.")
  }
  
  if (!is.numeric(burnin) || burnin < 0 || burnin != as.integer(burnin)) {
    stop("burnin must be a non-negative integer.")
  }
  
  if (!is.numeric(grid_size) || grid_size < 3 || grid_size != as.integer(grid_size)) {
    stop("grid_size must be an integer >= 3.")
  }
  
  if (!is.numeric(domain_length) || domain_length <= 0) {
    stop("domain_length must be a positive numeric value.")
  }
  
  if (!is.numeric(dt) || dt <= 0) {
    stop("dt must be a positive numeric value.")
  }
  
  if (!is.numeric(true_D) || true_D <= 0) {
    stop("true_D must be a positive numeric value.")
  }
  
  if (!is.numeric(true_rho) || true_rho <= 0) {
    stop("true_rho must be a positive numeric value.")
  }
  
  if (!is.logical(use_parallel)) {
    stop("use_parallel must be a logical value (TRUE or FALSE).")
  }
  
  if (!is.null(seed) && (!is.numeric(seed) || seed != as.integer(seed))) {
    stop("seed must be NULL or an integer value.")
  }
  
  # Check if reasonable
  if (burnin >= n_iter) {
    stop("burnin must be less than n_iter.")
  }
  
  if (noise_std > 0.5) {
    warning("High noise level (noise_std = ", noise_std, ") may lead to poor parameter identifiability.")
  }
  
  # Check if stable
  dx <- domain_length / (grid_size - 1)
  cfl <- true_D * dt / (dx * dx)
  if (cfl > 0.5) {
    warning(sprintf("CFL condition (%.3f) > 0.5 may cause numerical instability. Consider reducing dt.", cfl))
  }
  
  # Check if parallel works
  if (use_parallel && n_chains > 1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' not available. Falling back to sequential execution.")
      use_parallel <- FALSE
    }
  }
  
  # Execute function
  results <- .bayesian_tumor_growth_core(
    noise_std = noise_std,
    t_obs = t_obs,
    n_chains = n_chains,
    n_iter = n_iter,
    burnin = burnin,
    grid_size = grid_size,
    domain_length = domain_length,
    dt = dt,
    true_D = true_D,
    true_rho = true_rho,
    use_parallel = use_parallel,
    seed = seed
  )
  
  # Check results
  if (is.null(results)) {
    stop("Core function returned NULL results.")
  }
  
  required_components <- c("data", "chains", "D_mcmc", "rho_mcmc", 
                           "posterior_mean", "performance", "diagnostics")
  missing_components <- setdiff(required_components, names(results))
  if (length(missing_components) > 0) {
    stop("Core function missing required components: ", 
         paste(missing_components, collapse = ", "))
  }
  
  # Check convergence
  if (any(results$diagnostics$Rhat > 1.1, na.rm = TRUE)) {
    warning("Some Rhat values > 1.1 indicate potential convergence issues.")
  }
  
  # Check available samples 
  min_ess <- min(results$diagnostics$ESS)
  if (min_ess < 100) {
    warning(sprintf("Low effective sample size (min ESS = %.0f). Consider increasing n_iter.", min_ess))
  }
  
  return(results)
}

# Core function
.bayesian_tumor_growth_core <- function(noise_std, t_obs, n_chains, n_iter, burnin,
                                        grid_size, domain_length, dt, true_D, true_rho,
                                        use_parallel, seed) {
  # Generate synthesis data
  data <- .generate_tumor_data_internal(
    grid_size = grid_size,
    domain_length = domain_length,
    true_D = true_D,
    true_rho = true_rho,
    noise_std = noise_std,
    t_obs = t_obs,
    dt = dt,
    seed = seed
  )
  
  # Set prior
  prior <- list(
    D_mean_log = log(data$true_D),
    D_sd_log = 0.5,
    rho_mean_log = log(data$true_rho),
    rho_sd_log = 0.5
  )
  
  # Run multi-chain
  if (use_parallel && n_chains > 1) {
    cl <- parallel::makeCluster(min(n_chains, parallel::detectCores() - 1))
    parallel::clusterExport(cl, varlist = c("run_adaptive_mcmc_chain"), envir = environment())
    
    chains <- parallel::parLapply(cl, 1:n_chains, function(i) {
      run_adaptive_mcmc_chain(
        u_initial_r = data$u_initial,
        u_observed_r = data$u_observed,
        noise_std = data$noise_std,
        domain_length = data$domain_length,
        t_obs = data$t_obs,
        D_prior_mean_log = prior$D_mean_log,
        D_prior_sd_log = prior$D_sd_log,
        rho_prior_mean_log = prior$rho_mean_log,
        rho_prior_sd_log = prior$rho_sd_log,
        grid_size = length(data$u_initial),
        dt = data$dt,
        num_iterations = n_iter,
        burnin = burnin,
        seed = if (!is.null(seed)) seed + i else i * 1000
      )
    })
    
    parallel::stopCluster(cl)
  } else {
    chains <- vector("list", n_chains)
    for (i in seq_len(n_chains)) {
      chains[[i]] <- run_adaptive_mcmc_chain(
        u_initial_r = data$u_initial,
        u_observed_r = data$u_observed,
        noise_std = data$noise_std,
        domain_length = data$domain_length,
        t_obs = data$t_obs,
        D_prior_mean_log = prior$D_mean_log,
        D_prior_sd_log = prior$D_sd_log,
        rho_prior_mean_log = prior$rho_mean_log,
        rho_prior_sd_log = prior$rho_sd_log,
        grid_size = length(data$u_initial),
        dt = data$dt,
        num_iterations = n_iter,
        burnin = burnin,
        seed = if (!is.null(seed)) seed + i else i * 1000
      )
    }
  }
  
  # Trasnfer as coda objects
  D_mcmc <- coda::mcmc.list(lapply(chains, function(c) coda::mcmc(c$D)))
  rho_mcmc <- coda::mcmc.list(lapply(chains, function(c) coda::mcmc(c$rho)))
  
  # Diagonesis convergence
  Rhat_D <- tryCatch(coda::gelman.diag(D_mcmc)$mpsrf, error = function(e) NA)
  Rhat_rho <- tryCatch(coda::gelman.diag(rho_mcmc)$mpsrf, error = function(e) NA)
  ess_D <- coda::effectiveSize(D_mcmc)
  ess_rho <- coda::effectiveSize(rho_mcmc)
  
  # Perfermance 
  D_post_mean <- mean(sapply(chains, function(c) mean(c$D)))
  rho_post_mean <- mean(sapply(chains, function(c) mean(c$rho)))
  
  D_ae <- abs(D_post_mean - data$true_D)
  rho_ae <- abs(rho_post_mean - data$true_rho)
  D_re <- D_ae / data$true_D
  rho_re <- rho_ae / data$true_rho
  
  # Prediction
  u_pred <- simulate_pde_cpp(data$u_initial, D_post_mean, rho_post_mean,
                             data$domain_length, data$dt, data$t_obs,
                             length(data$u_initial))
  pred_mse <- mean((u_pred - data$u_true_final)^2)
  
  # CI coverage
  # Combine all chains and compute HPD interval
  D_all_chains <- unlist(lapply(D_mcmc, as.numeric))
  rho_all_chains <- unlist(lapply(rho_mcmc, as.numeric))
  
  D_CI <- coda::HPDinterval(coda::as.mcmc(D_all_chains))
  rho_CI <- coda::HPDinterval(coda::as.mcmc(rho_all_chains))
  
  # Extract lower and upper bounds (HPDinterval returns a matrix with "lower" and "upper" columns)
  cover_D <- as.numeric(D_CI[1, "lower"] <= data$true_D && data$true_D <= D_CI[1, "upper"])
  cover_rho <- as.numeric(rho_CI[1, "lower"] <= data$true_rho && data$true_rho <= rho_CI[1, "upper"])
  
  # Adaptive history
  adaptation_history <- lapply(chains, function(c) c$adaptation_history)
  
  list(
    data = data,
    chains = chains,
    D_mcmc = D_mcmc,
    rho_mcmc = rho_mcmc,
    posterior_mean = c(D = D_post_mean, rho = rho_post_mean),
    performance = list(
      absolute_error = c(D = D_ae, rho = rho_ae),
      relative_error = c(D = D_re, rho = rho_re),
      pred_mse = pred_mse,
      CI_coverage = c(D = cover_D, rho = cover_rho)
    ),
    diagnostics = list(
      Rhat = c(D = Rhat_D, rho = Rhat_rho),
      ESS = c(D = ess_D, rho = ess_rho),
      accept_rate = mean(sapply(chains, function(c) c$accept_rate))
    ),
    adaptation_history = adaptation_history
  )
}

# generate function with data
.generate_tumor_data_internal <- function(grid_size, domain_length, true_D, true_rho,
                                          noise_std, t_obs, dt, seed) {
  
  if (!is.null(seed)) set.seed(seed)
  
  x <- seq(0, domain_length, length.out = grid_size)
  u0 <- dnorm(x, mean = domain_length/2, sd = 3)
  u0 <- u0 / max(u0)  # normalization
  
  u_true <- simulate_pde_cpp(u0, true_D, true_rho, domain_length, dt, t_obs, grid_size)
  u_obs <- pmax(0, pmin(1, u_true + rnorm(length(u_true), 0, noise_std)))
    
  list(
    u_initial = u0,
    u_observed = u_obs,
    u_true_final = u_true,
    true_D = true_D,
    true_rho = true_rho,
    noise_std = noise_std,
    domain_length = domain_length,
    t_obs = t_obs,
    dt = dt,
    grid_size = grid_size
  )
}
