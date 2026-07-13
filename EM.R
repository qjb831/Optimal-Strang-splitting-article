library(MASS)
source('Brownian_motion.R',keep.source = TRUE)

# Simulerer repetitions af paths med EM-scheme
EM_sim <- function(F, G, x0, times, n_paths = 1, progress = TRUE) {
  d <- length(x0)
  N <- length(times) - 1
  dt <- diff(times)
  dW_list <- sim_BM(times, dim = d, M_reps = n_paths, get_increments = TRUE)
  
  paths <- vector("list", n_paths)
  
  # Initialize progress bar
  if (progress) {
    pb <- txtProgressBar(min = 0, max = n_paths, style = 3)
  }
  
  for (j in seq_len(n_paths)) {
    dW <- dW_list[[j]]
    X <- matrix(NA, N + 1, d)   # N+1 rows (include t0)
    X[1, ] <- x0
    for (i in 2:(N + 1)) {
      X[i, ] <- X[i - 1, ] + F(X[i - 1, ]) * dt[i - 1] + G(X[i - 1, ]) * as.numeric(dW[i - 1, ])
    }
    X <- data.frame(X)
    paths[[j]] <- X
    
    # update progress bar
    if (progress) setTxtProgressBar(pb, j)
  }
  
  if (progress) close(pb)
  
  return(paths)
}

EM_step <- function(F, G, x0, times, n=1000,progress=FALSE) {
  d <- length(x0)
  endpoints <- matrix(NA, n, d)  # store n endpoints
  
  for (i in 1:n) {
    path <- EM_sim_fast(F, G, x0, times, n_paths=1,progress=progress)[[1]]  # EM_sim returns a list
    endpoints[i, ] <- as.numeric(path[nrow(path), ]) # last row = endpoint
  }
  
  if (d==2) colnames(endpoints) <- c("x", "y")
  return(as.data.frame(endpoints))
}

EM_sim_fast <- function(F, G, x0, times, n_paths = 1, progress = FALSE) {
  d <- length(x0)
  N <- length(times) - 1
  dt <- diff(times)
  dW <- sim_BM(times, dim = d, M_reps = n_paths, get_increments = TRUE)
  
  X <- array(NA_real_, c(N + 1, d, n_paths))
  X[1, , ] <- x0
  
  if (progress) {
    pb <- txtProgressBar(min = 0, max = n_paths, style = 3)
    on.exit(close(pb), add = TRUE)  
  }
  
  # Outer loop over paths, inner loop over time
  for (j in seq_len(n_paths)) {
    for (i in 2:(N + 1)) {
      X[i, , j] <- X[i - 1, , j] +
        F(X[i - 1, , j]) * dt[i - 1] +
        G(X[i - 1, , j]) * as.numeric(dW[[j]][i - 1, ])
    }
    if (progress) setTxtProgressBar(pb, j)
  }
  
  new_names <- paste0("X", seq_len(d))
  
  df_list <- lapply(seq_len(n_paths), function(j) {
    df <- as.data.frame(X[, , j])
    colnames(df) <- new_names
    df
  })
  
  df_list
}


EM_step_fast <- function(F, G, x0, times, n = 1000) {
  sims <- EM_sim_fast(F, G, x0, times, n_paths = n)
  endpoints <- vapply(sims, function(path) as.numeric(path[nrow(path), ]), numeric(length(x0)))
  endpoints <- t(endpoints)
  colnames(endpoints) <- if (length(x0) == 2) c("x", "y") else paste0("x", seq_along(x0))
  as.data.frame(endpoints)
}

# Regner EM loglikelihood af en enkelt path
EM_loglik_path <- function(df, times, theta_drift, theta_diffusion, drift, diffusion) {
  n <- nrow(df)
  d <- ncol(df)

  dt <- diff(times)   # length = n-1
  drift_mat <- t(apply(df[-n, ], 1, function(x) drift(x, theta_drift)))
  
  sigma_vec <- diffusion(df[1, ], theta_diffusion)  
  vars <- matrix(rep(sigma_vec^2, each = n-1) * dt, nrow = n-1, ncol = d)
  vars <- pmax(vars, 1e-8)
  
  # observed increments
  X_diff <- df[-1, ] - df[-n, ]
  
  # log-likelihood contributions per step
  logliks <- -0.5*rowSums(log(vars))-0.5*rowSums((X_diff - drift_mat * dt)^2 / vars)
  
  sum(logliks)
}

# Regner samlet loglikelihood af liste af repetitions af paths 
EM_loglik_list <- function(list_of_df, times, theta_drift, theta_diffusion, drift, diffusion) {
  logliks <- vapply(
    list_of_df,
    function(df) {
      EM_loglik_path(df, times, theta_drift, theta_diffusion, drift, diffusion)
    },
    numeric(1)  # ensures numeric output
  )
  sum(logliks)
}

