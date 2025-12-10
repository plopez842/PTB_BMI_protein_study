# Purpose of Scripts----
# Prepare functions to run parallelization
compute_test_statistic_LRT <- function(model1,model2) {
  # Compute the deviance or any other appropriate test statistic
  
  A<-logLik(model2)
  B<-logLik(model1)
  
  testStat<-(-2)*(as.numeric(A)-as.numeric(B))
  # p.val<-pchisq(testStat,df=1,lower.tail=FALSE)
  return(testStat)
}



# Define a function to randomly permute the response variable
permute_response <- function(data, response_variable_name) {
  # Permute the response variable
  data[[response_variable_name]] <- sample(data[[response_variable_name]])
  return(data)
}


func1<-function(x,mat){ #x is the protein, matrix is comb_v0_long
  val<-mat[,x]
  bmi<-mat$bmi30
  type<-mat$Type
  
  temp_df<-as.data.frame(cbind(val,bmi,type))
  temp_df$val<-as.numeric(temp_df$val)
  
  
  form_v1<-as.formula(paste0("val ~ type + bmi + type*bmi"))
  gam_model_v1<-glm(form_v1,data=temp_df)
  
  #diff_AIC[diff_AIC$protein==i,]$AIC_1=AIC(gam_model_v1)
  
  form_v2<-as.formula(paste0("val ~ type + bmi"))
  
  gam_model_v2<-glm(form_v2,data=temp_df)
  
  
  
  observed_test_statistic <-compute_test_statistic_LRT(gam_model_v1,gam_model_v2)
  num_permutations <- 1000
  
  permuted_test_statistics <- numeric(num_permutations)
  for (i in 1:num_permutations) {
    # Permute the response variable in the data
    permuted_data <- temp_df
    permuted_data$bmi<-sample(permuted_data$bmi)
    # Fit a GAM model to the permuted data
    permuted_model_v1 <- glm(val ~ type + bmi + type*bmi, data=permuted_data)
    permuted_model_v2 <- glm(val ~ type + bmi,data=permuted_data)
    # Compute the test statistic for the permuted data
    permuted_test_statistics[i] <- compute_test_statistic_LRT( permuted_model_v1, permuted_model_v2)
  }
  p_value <- sum(permuted_test_statistics > observed_test_statistic) / num_permutations

  return(res)
}
