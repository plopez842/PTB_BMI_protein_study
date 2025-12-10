# Purpose of Script ---- 
# Figures for Figure 2C

# Load Variables----
set.seed(1010)
source("../helpful_functions/Systems_Serology_PAL.R")
source("../helpful_functions/fig2d_univ_plot.R")

suppressPackageWarnings({
  library(ggpubr)
  library(dplyr)
})

# Setup of Data Frames----

## Prepare DF ----
patient_info_levels = readRDS("../output_RDS_objs/merged_metadata_and_protein_V1.rds")
proteins<-colnames(patient_info_levels)[grepl("seq.",colnames(patient_info_levels))]
merged_data_mod<-patient_info_levels

## Fix some of the variable names ---
merged_data_mod$Type<-gsub("meta_sptb","sptb",merged_data_mod$Type)
merged_data_mod$Type<-gsub("mtpb","mptb",merged_data_mod$Type)

merged_data_mod$bmi30<-gsub("<30","less30",merged_data_mod$bmi30)
merged_data_mod$bmi30<-gsub(">=30","greater30",merged_data_mod$bmi30)


## Prepare the X for PLSDA----
rownames(merged_data_mod)<-paste0("studyid_",merged_data_mod$studyid)
X_unfilt_data<-merged_data_mod[,grepl("seq.",colnames(merged_data_mod))]
X_scale<-scale(X_unfilt_data,center=TRUE,scale=TRUE)

## Prepare the Y for PLSDA----
y_var<-merged_data_mod$Type
y_var<-factor(y_var,levels=c("ctrl","mptb","sptb"))


# Setup Model parameters----
opts = list(n_folds = 5, pt_trials =100, rf_trials=100,save = TRUE, compare_pred = "y")
vals_trial<-list()


method = list(select=feature_select_iter,
              train = train_ropls,
              predict = predict_ropls,
              score = score_accuracy)
# Run Validation----

vals <- cross_validation_unpaired(X_scale,y_var, method, opts,n_trials=10)


# Plot results----
plt_validation<- visualize_validate(vals, options = list(y_label = "accuracy")) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 12),
        axis.text.x=element_text(size=18))
# Save Results ----
ggsave(filename = "Fig2C_model_val.pdf",plt_val_TRIAL1,units="in",
       path="../final_plots/figure2",
       width=4,
       height=4)

# Save CSV----
all_comb = data.frame(Round = 1,
                      CV = rep(vals[[1]]$cv_score,100),
                      random_feat = vals[[1]]$rf_scores,
                      permuted_labels = vals[[1]]$pt_scores)

for(i in 2:n_trials){
  temp_df =  data.frame(Round = i,
                        CV = rep(vals[[i]]$cv_score,100),
                        random_feat = vals[[i]]$rf_scores,
                        permuted_labels = vals[[i]]$pt_scores)
  
  all_comb = rbind(all_comb,temp_df)
  
}
write.csv(all_comb,"../Supplemental_data_revisions/fig2c_REAL_model_val.csv")
