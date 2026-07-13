sim_BM <- function(times, dim = 1, M_reps = 1,get_increments = FALSE) {
  N <- length(times) - 1       # 
  dt <- diff(times)     
  
  reps <- vector("list", M_reps)
  for (i in seq_len(M_reps)) {
    sd_vec <- sqrt(rep(dt, each = dim))
    
    dW <- rnorm(dim*N, mean = 0, sd = sd_vec)
    dW <- matrix(dW, nrow = dim, ncol = N) 
    if (get_increments == TRUE) {
      dW  <- as.data.frame(t(dW))
      colnames(dW) <- gsub("V", "dW", colnames(dW))
      reps[[i]] <- dW
      }
    else{
      B <- apply(dW, 1, cumsum)  
      
      if (is.vector(B)) {
        # handle dim == 1 where apply returns a vector
        B <- matrix(B, nrow = dim, ncol = N)
      }
      B <- rbind(rep(0, dim),B) 
      B <- as.data.frame(B)
      colnames(B) <- gsub("V", "W", colnames(B))
      reps[[i]] <- B
    }
  }
  return(reps)
}

