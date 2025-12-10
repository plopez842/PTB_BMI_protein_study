# Purpose of Script----
# Run graph lasso on mPTB vs Corr matrix

# Load Functions----

functions_path_dir="../helpful_functions/fig4_PCN"
r_files <- list.files(path = functions_path_dir, pattern = "\\.R$", full.names = TRUE)
sapply(r_files,source)


# Load Objcts----
df_int <- readRDS("../output_RDS_objs/fig5_mPTB_for_PCN.rds")

# Graph Lasso Parameters----

#variable input
daterun = "20250304"
min = 0.02
max=0.2
interval = 0.01

#lambda thresholds
cutoff_min = 0.9
cutoff_max = 0.99
cutoff_interval = 0.01
#how many perm tests to do
perm_count = 100

#save the user parameters
user_param = list()
user_param[['lambda_min']] = min
user_param[['lambda_max']] = max
user_param[['lambda_interval']] = interval
user_param[['cutoff_min']] = cutoff_min
user_param[['cutoff_max']] = cutoff_max
user_param[['cutoff_interval']] = cutoff_interval
user_param[['perm_count']] = perm_count

# Run Graph Lasso----
test_rows<-colnames(df_int)
data_matrix=df_int
data_matrix = preprocess_data(data_matrix)
lambda_range=seq(min,max,by=interval)
glasso_output<-stable_gLASSO(data_matrix,test_rows,
                             lambda_list=lambda_range,
                             scale=FALSE) #this parameters scale internally)

samplecov<-cov(data_matrix[, test_rows], use = "complete.obs") # not scaled 

full_stabmat<-glasso_output$stability_matrix
colnames(full_stabmat)<-lambda_range

# Adjust----

stabmat=full_stabmat[,as.character(seq(min,max,by=interval))]

p=length(test_rows)
adj_mat_list=list()

#this is to setup varying cutoff ranges
num_iter=perm_count#this is how many iterations to test
cutoff_range<-seq(cutoff_min,cutoff_max,by=cutoff_interval)

comp_mat=gen_comp_mat(num_iter = num_iter,
                      cutoff_range=cutoff_range)
# add in max and mean real-est error
comp_mat$mean_error = NA
comp_mat$max_error = NA
for(j in cutoff_range){
  adjmat<-matrix(as.numeric(apply(stabmat, 1, max) > j), p, p)
  #determine edge density
  edge_density=mean(adjmat[upper.tri(adjmat)])
  comp_mat[comp_mat$cutoff==j,]$edge_density=edge_density
  
  #compare the correlation matrix only the unadjusted values
  
  ##fix the real
  samplecov_sub<-sub_cov_adj_mat(adjmat,samplecov)
  
  ##get the estimated cov
  est_cov_real<-stable_partcorr(adjmat,samplecov)
  ###get the cov values of the indirect correlates
  sub_est_cov_real<-sub_cov_adj_mat(adjmat,est_cov_real$covarmat)
  
  
  #setup r2 
  ## this is to determien the position of the INDIRECT CORRELATES, starting with the real cov mat
  adj_mat_zeros = get_off_diagonal(adjmat)
  adj_mat_zeros = adj_mat_zeros==0
  
  ##real
  diag_sample <- get_off_diagonal(samplecov_sub)
  diag_sample <- diag_sample[adj_mat_zeros]
  
  ##est
  diag_est <- get_off_diagonal(sub_est_cov_real)
  diag_est <- diag_est[adj_mat_zeros]
  
  r2_val <- summary(lm(diag_sample ~ diag_est))$r.squared
  
  comp_mat[comp_mat$var=="act" & comp_mat$cutoff==j,]$R2=r2_val
  mean_error = mean(abs(diag_sample-diag_est))
  max_error = max(abs(diag_sample-diag_est))
  
  comp_mat[comp_mat$var=="act" & comp_mat$cutoff==j,]$mean_error = mean_error
  comp_mat[comp_mat$var=="act" & comp_mat$cutoff==j,]$max_error = max_error
  
  #setup for perm
  perm_comp_mat = comp_mat[comp_mat$cutoff==j & comp_mat$var!="act",]
  for(i in 1:num_iter){
    print(i)
    rand_adjmat_shuffle<-shuffle_symmetric(adjmat)
    
    perm_par_corr<-stable_partcorr(rand_adjmat_shuffle,samplecov)
    perm_cov_mat=perm_par_corr$covarmat
    
    sub_est_cov_shuffle<-sub_cov_adj_mat(adjmat,perm_cov_mat) #only look at estimated partial correlates
    
    diag_perm <- get_off_diagonal(sub_est_cov_shuffle)
    diag_perm <- diag_perm[adj_mat_zeros]
    
    r2_val_perm=summary(lm(diag_sample ~ diag_perm))$r.squared
    mean_error_perm = mean(abs(diag_sample-diag_perm))
    max_error_perm = max(abs(diag_sample-diag_perm))
    perm_comp_mat[perm_comp_mat$perm_RD==paste0("perm_RD_",i),]$R2 = r2_val_perm
    perm_comp_mat[perm_comp_mat$perm_RD==paste0("perm_RD_",i),]$mean_error =  mean_error_perm
    perm_comp_mat[perm_comp_mat$perm_RD==paste0("perm_RD_",i),]$max_error =  max_error_perm
  }
  comp_mat[comp_mat$cutoff==j & comp_mat$var!="act",] = perm_comp_mat
}

# Save RDS objects----

#save comp matrix
filename_r2_act_perm = paste0(daterun,"_mPTB_r2_real_perm_lambda_min_",min,"_lambda_max_",max,".Rds")
saveRDS(comp_mat,paste0("../output_RDS_objs/Fig5_PCN_int_objects/",filename_r2_act_perm))

#save glasso results
filename_stab_mat = paste0(daterun,"_mPTB_glasso_stab_mat_lambda_min_",min,"_lambda_max_",max,".Rds")
saveRDS(full_stabmat,paste0("../output_RDS_objs/Fig5_PCN_int_objects/",filename_stab_mat))

#save user results 
filename_user_param = paste0(daterun,"mPTB_USER_PARAM_lambda_min_",min,"_lambda_max_",max,".RData")
save(user_param,file=paste0("../output_RDS_objs/Fig5_PCN_int_objects/",filename_user_param))

