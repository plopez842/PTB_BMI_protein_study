# Purpose of Script----
# Prepare variables for the PCN


## Proteins included in the PCN----
# - 40 proteins influenced solely based on BMI (BMI SPECIFIC) - use pairwise comparisons
# - 40 proteins based on preterm birth (PTB SPECIFIC) - use LV1 loadings of PLSDA model)
# - 40 proteins from LRT (mPTB vs Ctrl)
# - 40 proteins from LRT (msTB vs Ctrl)
# Note: Some aptamers were duplicates (targeting the same protein. Chose to continue with the highest F_statistic
# Note: Also some were hits in multiple criterea

# Load Libraries----
suppressPackageStartupMessages({
  library(dplyr)
})

# Load Variables----

patient_info_levels = readRDS("../output_RDS_objs/merged_metadata_and_protein_V1.rds")
proteins<-colnames(patient_info_levels)[grepl("seq.",colnames(patient_info_levels))]
merged_data_mod<-patient_info_levels

## load Patient and fix columns----
merged_data_mod$Type<-gsub("meta_sptb","sptb",merged_data_mod$Type)
merged_data_mod$Type<-gsub("mtpb","mptb",merged_data_mod$Type)

merged_data_mod$bmi30<-gsub("<30","less30",merged_data_mod$bmi30)
merged_data_mod$bmi30<-gsub(">=30","greater30",merged_data_mod$bmi30)
merged_data_mod$group<-paste0(merged_data_mod$Type,"_",merged_data_mod$bmi30)
group<-merged_data_mod$group

## Work on Apt DF ----
protein_apt_info = readRDS("../output_RDS_objs/Somalogic_protein_apt_info.rds")
max_entrez<-protein_apt_info
max_entrez<-max_entrez %>%
  group_by(Single_Entrez) %>%
  slice_max(f_stat) %>%
  ungroup()
protein_apt_info$max_f_stat_single_entrez="N"
protein_apt_info[protein_apt_info$AptName %in% max_entrez$AptName,]$max_f_stat_single_entrez="Y"

## Load the Protein Data from other measurements----

### BMI-specific signal----
pairwise_comp_ALL <- readRDS(file="../output_RDS_objs/pairwise_comp_volc_plot_list_df.rds")
pairwise_comp<-pairwise_comp_ALL$g30Vsl30

sub_pairwise_comp<-subset(pairwise_comp,adj.P.Val <= 0.05)
sub_pairwise_comp<-sub_pairwise_comp[order(abs(sub_pairwise_comp$logFC),decreasing=TRUE),]
sub_pairwise_comp<-sub_pairwise_comp[1:40,]

###PTB-specific signal----
PTB_Ctrl_corr_vals<-readRDS(file="../output_RDS_objs/fig2_PLSDA_mPTB_sPTB_Ctrl.rds")
PTB_Ctrl_corr_vals <- PTB_Ctrl_corr_vals[order(abs(PTB_Ctrl_corr_vals$rho),decreasing=TRUE),]
PTB_Ctrl_corr_vals<-PTB_Ctrl_corr_vals[1:40,]

### LRT-mPTB----
LRT_mPTB<-readRDS(file="../output_RDS_objs/fig4c_LRT_results_mPTB_Ctrl.rds")
LRT_mPTB<-subset(LRT_mPTB,pval_LRT_actual < 0.005)
LRT_mPTB<-LRT_mPTB[order(abs(LRT_mPTB$pval_LRT_actual),decreasing=FALSE),]
LRT_mPTB<-LRT_mPTB[1:40,]

### LRT-sPTB----
LRT_sPTB<-readRDS(file="../output_RDS_objs/fig4d_LRT_results_sPTB_Ctrl.rds")
LRT_sPTB<-subset(LRT_sPTB,pval_LRT_actual < 0.005)
LRT_sPTB<-LRT_sPTB[order(abs(LRT_sPTB$pval_LRT_actual),decreasing=FALSE),]
LRT_sPTB<-LRT_sPTB[1:40,]

# Prep Proteins of Int----
PROTEINS_INT <- unique(rownames(sub_pairwise_comp),
                  PTB_Ctrl_corr_vals$protein,
                  LRT_mPTB$protein,
                  LRT_sPTB$protein)

## Handle Duplciates----
sub_protein_apt_info<-subset(protein_apt_info,AptName %in% PROTEINS_INT)
filtered_protein_apt_info <- protein_apt_info %>%
  group_by(EntrezGeneSymbol) %>%          # Group by the EntrezGeneSymbol column
  slice_max(f_stat, n = 1, with_ties = FALSE) %>%  # Keep the row with the highest f_stat
  ungroup()   


PROTEINS_INT <- filtered_protein_apt_info$AptName

# Prep DF of samples
## mPTB
PTB_int = "mptb"
sub_type<-subset(merged_data_mod,Type %in% c("ctrl",PTB_int))
sub_type$Type2<-factor(sub_type$Type,levels=c("ctrl",PTB_int))
sub_type$Type2<-factor(sub_type$Type2,labels=c(0,1))
sub_type$Type2<-as.numeric(as.vector(sub_type$Type2))
sub_mptb = sub_type

## sPTB
PTB_int = "sptb"
sub_type<-subset(merged_data_mod,Type %in% c("ctrl",PTB_int))
sub_type$Type2<-factor(sub_type$Type,levels=c("ctrl",PTB_int))
sub_type$Type2<-factor(sub_type$Type2,labels=c(0,1))
sub_type$Type2<-as.numeric(as.vector(sub_type$Type2))
sub_sptb = sub_type

## subset df to desired columns
test_rows = c("bmi1","Type2",PROTEINS_INT)

# Prepare protein info df----


proteins_int_info<-data.frame(apt_name = filtered_protein_apt_info$AptName,
                              node_name = filtered_protein_apt_info$Single_gene,
                              protein_or_covariate = "protein",
                              test1_BMI_only=0,
                              test2_PTB_Ctrl_PLSDA=0,
                              test3_LRT_mPTB=0,
                              test4_LRT_sPTB=0,
                              test1_val=NA,
                              test2_val=NA,
                              test3_val=NA,
                              test4_val=NA)

list_tests=list(sub_pairwise_comp,
                PTB_Ctrl_corr_vals,
                LRT_mPTB,
                LRT_PTB)
names(list_tests)<-c("test1","test2","test3","test4")

for(i in names(list_tests)){
  df_int = list_tests[[i]]
  
  
  #first acquire proteins int and values
  if(i == "test1"){
    temp_proteins_int = rownames(df_int)
    temp_vals = df_int$logFC
  }else{
    temp_proteins_int = df_int$protein
    if(i=="test2"){
      temp_vals = df_int$rho
    }else{
      temp_vals = df_int$log2_diff #not actually used to threshold but interesting to see 
    }
  }
  temp_df<-data.frame(protein_int = temp_proteins_int,
                      val_int = temp_vals)
  
  
  #now working on inserting
  if(i=="test1"){
    proteins_int_info[proteins_int_info$apt_name %in% temp_df$protein_int,]$test1_BMI_only=1
    proteins_int_info$test1_val<-pairwise_comp[match(proteins_int_info$apt_name,rownames(pairwise_comp)),]$logFC
    
  }else if(i=="test2"){
    proteins_int_info[proteins_int_info$apt_name %in% temp_df$protein_int,]$test2_PTB_Ctrl_PLSDA=1
    proteins_int_info$test2_val <-temp_df$val_int[match(proteins_int_info$apt_name,temp_df$protein_int)]
  }else if(i=="test3"){
    proteins_int_info[proteins_int_info$apt_name %in% temp_df$protein_int,]$test3_LRT_mPTB=1
    proteins_int_info$test3_val <-temp_df$val_int[match(proteins_int_info$apt_name,temp_df$protein_int)]
  }else{
    proteins_int_info[proteins_int_info$apt_name %in% temp_df$protein_int,]$test4_LRT_sPTB=1
    proteins_int_info$test4_val <-temp_df$val_int[match(proteins_int_info$apt_name,temp_df$protein_int)]
  }
}


df_proteins_info<-data.frame(protein = proteins_int_info$node_name,
                             Significant_in_which_test= NA)

for(i in 1:nrow(proteins_int_info)){
  temp_row = proteins_int_info[i,]
  
  blank_name = c()
  if(temp_row$test1_BMI_only==1){
    blank_name = c(blank_name,'Sig Between BMI groups (<30 vs >30 kg/m2)')
  }
  
  if(temp_row$test2_PTB_Ctrl_PLSDA==1){
    blank_name = c(blank_name,'Sig Between PTB vs Ctrl - LV1 of 3-way PLSDA')
  }
  
  if(temp_row$test3_LRT_mPTB){
    blank_name = c(blank_name,'Sig in LRT test - mPTB vs Ctrl')
  }
  
  if(temp_row$test4_LRT_sPTB){
    blank_name = c(blank_name,'Sig in LRT test - sPTB vs Ctrl')
  }
  
  merged_name = paste(blank_name,collapse = "; ")
  
  df_proteins_info[i,]$Significant_in_which_test = merged_name
}


## save RDS
saveRDS(sub_mptb[,test_rows],file = "../output_RDS_objs/fig5_mPTB_for_PCN.rds")
saveRDS(sub_sptb[,test_rows],file = "../output_RDS_objs/fig5_sPTB_for_PCN.rds")
saveRDS(df_proteins_info,file = "../output_RDS_objs/fig5_protein_info.rds")
