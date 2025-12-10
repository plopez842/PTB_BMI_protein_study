stable_partcorr <- function(adjmat, samplecov, tol = 1e-6) {
  n <- nrow(adjmat)
  diag(adjmat) <- 0  # Ensure no self-loops in adjacency matrix
  
  precmat <- matrix(0, nrow = n, ncol = n)
  W <- samplecov
  W_iter <- samplecov  # Initialize with sample covariance matrix
  converg_flag <- FALSE
  count <- 0
  
  while (!converg_flag) {
    count <- count + 1
    for (j in 1:n) {
      idx <- rep(TRUE, n)  # Select all variables
      idx[j] <- FALSE  # Exclude j from neighbors
      
      W_11 <- W_iter[idx, idx, drop = FALSE]
      beta <- rep(0, sum(idx))  # Initialize beta of appropriate size
      neighb_id <- adjmat[idx, j] != 0  # Boolean indicating neighbors
      
      neighb_idx <- which(neighb_id)  # Indices of neighbors
      W_11s <- W_iter[idx & adjmat[, j] != 0, idx & adjmat[, j] != 0, drop = FALSE]
      s_12s <- samplecov[idx & adjmat[, j] != 0, j]
      
      if (nrow(W_11s) > 0) {
        # Solve for beta coefficients
        betas <- try(solve(W_11s, s_12s), silent = TRUE)
        
        # Handle the case where solve fails due to singularity
        if (inherits(betas, "try-error")) {
          stop("Matrix is singular or near-singular. No regularization is used as per request.")
        }
        
        beta[neighb_id] <- betas  # Fill in the non-zero elements of beta
      }
      
      # Update W_iter with the new estimates
      W_iter[j, idx] <- as.vector(W_11 %*% beta)
      W_iter[idx, j] <- as.vector(W_11 %*% beta)
    }
    
    # Check for convergence
    
    #print(max(abs(W - W_iter)))
    if (max(abs(W - W_iter)) < tol) {
      covarmat <- W_iter
      precmat <- try(solve(covarmat), silent = TRUE)  # Inverse of the covariance matrix
      
      # Handle case where covariance is singular
      if (inherits(precmat, "try-error")) {
        stop("Covariance matrix is singular and cannot be inverted.")
      }
      
      converg_flag <- TRUE
    } else {
      W <- W_iter  # Update W for the next iteration
    }
  }
  
  cat('Number of iterations:', count, '\n')
  return(list(precmat=precmat,covarmat=covarmat,count=count))
}
# 


