suppressPackageStartupMessages({
    library(glasso) # For graphical LASSO
    library(foreach) # For parallel looping (optional)
    library(doSNOW) # For the snow parallel backend (optional)
    library(RhpcBLASctl) # For BLAS threading control
    library(MASS) # For multivariate normal functions
   # library(emdbook)
})

stable_gLASSO <- function(DataFull, varid, lambda_list, k=5,
                          scale=FALSE,
                          n=100,
                          parallelize = FALSE, 
                          num_cores = 2) {
  # DataFull: A data frame containing the data to be modeled
  # varid: A vector indicating the columns of DataFull to use for analysis
  # lambda_list: A list containing the LASSO lambda parameters
  # k: k-fold cross-validation parameter
  # n: Number of data subsamples to take
  # parallelize: Boolean to determine whether to run the code in parallel (default = FALSE)
  # num_cores: Number of cores to use for parallel processing if parallelize = TRUE (default = 2)
  
  # Set BLAS threads to 1 to avoid conflicts with parallelization
  blas_set_num_threads(1)
  
  # Preprocess the data
  
  ## do we need to scale?
  if(scale){
    data <- preprocess_data(DataFull[, varid])
  }
  else{
    data <- (DataFull[, varid])
    }
  
  p <- ncol(data)
  lambda_res <- length(lambda_list)
  zero_cutoff <- 1e-4
  
  # Initialize matrices
  LLH_matrix <- matrix(0, nrow = n, ncol = lambda_res)
  stability_matrix <- matrix(0, nrow = p^2, ncol = lambda_res)
  
  # If parallelization is enabled, initialize the cluster
  if (parallelize) {
    cl <- makeCluster(num_cores, type = "SOCK")
    registerDoSNOW(cl)
  }
  
  # Iterate over lambda values
  for (i in 1:lambda_res) {
    stability_counter <- matrix(0, nrow = p^2, ncol = n)
    
    if (parallelize) {
      # Use foreach for parallel subsampling
      results <- foreach(j = 1:n, .combine = 'rbind', .packages = c("glasso", "MASS")) %dopar% {
        split_result <- splitdata(data, k)
        testdata <- split_result$testdata
        traindata <- split_result$traindata
        
        # Run graphical LASSO on the training set using glasso with the specified lambda
        glasso_result <- glasso(cov(traindata), rho = lambda_list[i])
        Theta_train <- glasso_result$wi  # The inverse covariance matrix (precision matrix)
        
        # Stability counter: whether the edge is non-zero
        stability_vec <- abs(as.vector(Theta_train)) > zero_cutoff
        
        # Evaluate the log-likelihood on the test set
        W_train <- cov2cor(cov(traindata))  # Convert to correlation matrix for LLH
        LLH_value <- -gaussianLLH(testdata, W_train)
        
        # Return results: stability vector and log-likelihood value
        list(stability_vec = stability_vec, LLH_value = LLH_value)
      }
      
      # Extract and combine results from foreach
      stability_counter <- do.call(cbind, lapply(results, function(x) x$stability_vec))
      LLH_matrix[, i] <- unlist(lapply(results, function(x) x$LLH_value))
      
    } else {
      # Sequential subsampling if not parallelized
      for (j in 1:n) {
        split_result <- splitdata(data, k)
        testdata <- split_result$testdata
        traindata <- split_result$traindata
        
        # Run graphical LASSO on the training set using glasso with the specified lambda
        glasso_result <- glasso(cov(traindata), rho = lambda_list[i])
        Theta_train <- glasso_result$wi  # The inverse covariance matrix (precision matrix)
        
        # Stability counter: whether the edge is non-zero
        stability_counter[, j] <- abs(as.vector(Theta_train)) > zero_cutoff
        
        # Evaluate the log-likelihood on the test set
        W_train<-glasso_result$w #estimated covariance matrix for LLH 
        LLH_matrix[j, i] <- -gaussianLLH(testdata, W_train)
      }
    }
    
    # Calculate stability matrix
    stability_matrix[, i] <- rowMeans(stability_counter)
  }
  
  # If parallelization was used, stop the cluster
  if (parallelize) {
    stopCluster(cl)
  }
  
  # Reset BLAS threads to default
  blas_set_num_threads(RhpcBLASctl::blas_get_num_procs())
  
  # Return the results
  list(stability_matrix = stability_matrix, LLH_matrix = LLH_matrix)
}

# Preprocessing function
preprocess_data <- function(data_in) {
  null_id <- apply(is.na(data_in), 1, any)
  data <- scale(data_in[!null_id, ])
  return(data)
}

# Split data function for k-fold cross-validation
splitdata <- function(data, k) {
  n_sample <- nrow(data)
  test_id <- sample(1:n_sample, floor(n_sample / k))
  train_id <- setdiff(1:n_sample, test_id)
  list(testdata = data[test_id, ], traindata = data[train_id, ])
}

# Gaussian log-likelihood function
gaussianLLH <- function(sampledata, cov_train) {
  # sampledata should have zero mean
  #llh_values <- dmvnorm(sampledata, mu = rep(0, ncol(sampledata)), Sigma = covarmat, log = TRUE)
  #sampledata: Data from test set, assumed to have zero mean
  # cov_train: Regularized covariance matrix estimated by graphical lasso
  
  n <- nrow(sampledata)
  p <- ncol(sampledata)
  
  # Compute the log determinant of the covariance matrix
  log_det_cov <- determinant(cov_train, logarithm = TRUE)$modulus
  
  # Compute the quadratic term: (X^T Σ^(-1) X) for the multivariate Gaussian
  quad_term <- sum(diag(sampledata %*% solve(cov_train) %*% t(sampledata)))
  
  # Log-likelihood for multivariate Gaussian distribution
  LLH <- -0.5 * (n * (p * log(2 * pi) + log_det_cov) + quad_term)
  return(LLH)
}


get_PC_glasso_matrix<-function(adjmat,cov_mat){
  red_cov_mat = cov_mat
  for(i in 1:nrow(adjmat)){
    for(j in 1:ncol(adjmat)){
      if(adjmat[i,j]==0){
        red_cov_mat[i,j]=0
      }
    }
  }
  return(red_cov_mat)
}
