
gen_comp_mat<-function(num_iter,cutoff_range){
  cutoff_col=rep(cutoff_range,(num_iter)+1) #column1, adding 1 to account for actual val set
  edge_col<-cutoff_col #column 2
  
  act_var=rep('act',length(cutoff_range))
  #perm=rep('perm',length(cutoff_range))
  perm_var=rep("perm",length(cutoff_range)*num_iter)
  variable=c(act_var,perm_var) #column3
  
  perm_rd=paste0(perm_var,"_RD_",rep(1:num_iter,each=length(cutoff_range)))
  RD_col=c(act_var,perm_rd) #column 4
  
  r2_val=rep(NA,length(cutoff_col)) #column5
  
  
  comp_mat=data.frame(cutoff=cutoff_col,
                      edge_density=edge_col,
                      var=variable,
                      perm_RD=RD_col,
                      R2=r2_val)
  return(comp_mat)}




comp_cutoff_ranges<-function(stabmat,cutoff_range,comp_mat,samplecov){
  p=sqrt(nrow(stabmat))
  for(j in cutoff_range){
    #work with estimated cov
    adjmat<-matrix(as.numeric(apply(stabmat, 1, max) > j), p, p)
    est_cov_mat <- stable_partcorr(adjmat,samplecov)
    sub_est_cov_mat <- sub_cov_adj_mat(adjmat,est_cov_mat$covarmat)
    
    #get real
    samplecov_sub<-sub_cov_adj_mat(adjmat,samplecov)
    r2_est_act_cov(sub_est_cov_mat,samplecov_sub)
    #adj_mat_list[[paste0("cutoff_",j)]]=adjmat
    #determine edge density
    edge_density=mean(adjmat[upper.tri(adjmat)])
    comp_mat[comp_mat$cutoff==j,]$edge_density=edge_density
    
    
    r2_val=r2_adj_mat(sub_est_cov_mat,samplecov_sub)
    #print(r2_val)
    
    comp_mat[comp_mat$var=="act" & comp_mat$cutoff==j,]$R2=r2_val
  }
  return(comp_mat)
  
}


extract_numbers <- function(text) {
  as.numeric(unlist(regmatches(text, gregexpr("[0-9]+\\.[0-9]+", text))))
}

#this is for when generating 

shuffle_symmetric <- function(mat) {
  n <- nrow(mat)
  
  # Shuffle only the upper triangular part excluding the diagonal
  upper_tri_indices <- which(upper.tri(mat))
  shuffled_vals <- sample(mat[upper_tri_indices])
  
  # Assign the shuffled values to the upper and lower triangle
  mat[upper_tri_indices] <- shuffled_vals
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  
  # Ensure diagonal remains 1
  diag(mat) <- 1
  
  return(mat)
}

get_off_diagonal <- function(mat) {
  return(mat[upper.tri(mat, diag = FALSE)])
}

r2_est_act_cov<-function(est_cov,act_cov){

  actual_off_diag <- get_off_diagonal(act_cov)
  estimated_off_diag <- get_off_diagonal(est_cov)
  r_squared <- summary(lm(estimated_off_diag ~ actual_off_diag))$r.squared
  
  return(r_squared)
}

# the function below makes the partial correlates 0 if they have a 1 in adj mat
sub_cov_adj_mat<-function(adjmat,cov_mat){
  red_cov_mat = cov_mat
  for(i in 1:nrow(adjmat)){
    for(j in 1:ncol(adjmat)){
      if(adjmat[i,j]==1){
        red_cov_mat[i,j]=0
      }
    }
  }
  return(red_cov_mat)
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

