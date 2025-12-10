# Purpose of Script----
## Analysis and results of Pertubation Analysis

# Load Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Load Variables

## Load Precision matrix
mPTB_precision_matrix = readRDS(file = "../output_RDS_objs/Fig5_PCN_int_objects/Fig5_mPTB_precision_matrix.RDS")
sPTB_precision_matrix = readRDS(file = "../output_RDS_objs/Fig5_PCN_int_objects/Fig5_sPTB_precision_matrix.RDS")
## Work on Apt DF 
protein_apt_info = readRDS("../output_RDS_objs/Somalogic_protein_apt_info.rds")
max_entrez<-protein_apt_info
max_entrez<-max_entrez %>%
  group_by(Single_Entrez) %>%
  slice_max(f_stat) %>%
  ungroup()
protein_apt_info$max_f_stat_single_entrez="N"
protein_apt_info[protein_apt_info$AptName %in% max_entrez$AptName,]$max_f_stat_single_entrez="Y"

# Apply Perturbation just for mPTB----
prec_mat = mPTB_precision_matrix 
nodes_in_prec_matrix = colnames(prec_mat) #this also incldues the bmi and Type2 cols, important for indexing
for(i in feat_int){
  sub_prec_matrix<-prec_mat[nodes_in_prec_matrix !=i,nodes_in_prec_matrix !=i]
  new_corr_mat<-solve(sub_prec_matrix)
  new_pred_bmi_corr<-new_corr_mat['bmi1','Type2']
  predicted_corr_bmi_ptb[predicted_corr_bmi_ptb$protein==i,]$est_corr_val=new_pred_bmi_corr
}


## Merge with Apt Info----
sub_protein_apt_info<-subset(protein_apt_info,AptName %in% feat_int)
gene_int<-sub_protein_apt_info[match(feat_int,protein_apt_info$AptName),]$EntrezGeneSymbol

predicted_corr_bmi_ptb$Gene_Name = gene_int

## Take Diff Between real vs actual
predicted_corr_bmi_ptb$difference =  est_cov_bmi_ptb - (predicted_corr_bmi_ptb$est_corr_val)

# Plot (mPTB)----

## Cutoffs determined from Permutation Analysis done by Luka 
FDR_cutoff_2 = -0.06719684114495852
FDR_cutoff_1 =0.08664995053599511
high_diff_proteins<-subset(predicted_corr_bmi_ptb,difference > FDR_cutoff_1)$Gene_Name
others_to_keep = c("ENOX2","CRP","USP11","EVA1C","LAS2","GMEB2") # Results from sPTB (below)

all_together = c(high_diff_proteins,others_to_keep)
# Subset your data first
plot_data <- subset(predicted_corr_bmi_ptb, Gene_Name %in% all_together)

# Reverse sort factor levels alphabetically
plot_data$Gene_Name <- factor(plot_data$Gene_Name, levels = sort(unique(plot_data$Gene_Name), decreasing = TRUE))

plot_data_mPTB = plot_data 
bar_plot_pert_mPTB<-ggplot(plot_data,aes(x=Gene_Name,y=difference))+
  geom_col(stat="identity",width = 0.5,color = "black",fill= "gray")+
  coord_flip()+
  #geom_text(aes(label = gene_label),hjust=-0.2,size=2,vjust=0.1)+
  theme_classic()+
  ggtitle("BMI-mPTB diff. corr.\npost-perturb")+
  #adsd lines
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = FDR_cutoff_1,linetype = "dashed",color = "red")+
  geom_hline(yintercept = FDR_cutoff_2,linetype = "dashed",color = "red")+
  xlab("Perturbed Protein")+
  ylab("Diff. of Unperturb - Perturb Corr\nBMI to PTB")+
  scale_y_continuous(
    breaks = seq(-0.15,0.2,by=0.1),
    limits = c(-0.2,0.2)
  )+
  theme(#axis.text.y = element_blank(),
    plot.title = element_text(size=16,hjust=0.5,face="bold"),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.title = element_text(face="bold",size=12))


# Apply Perturbation just for mPTB----
prec_mat = mPTB_precision_matrix 
nodes_in_prec_matrix = colnames(prec_mat) #this also incldues the bmi and Type2 cols, important for indexing
for(i in feat_int){
  sub_prec_matrix<-prec_mat[nodes_in_prec_matrix !=i,nodes_in_prec_matrix !=i]
  new_corr_mat<-solve(sub_prec_matrix)
  new_pred_bmi_corr<-new_corr_mat['bmi1','Type2']
  predicted_corr_bmi_ptb[predicted_corr_bmi_ptb$protein==i,]$est_corr_val=new_pred_bmi_corr
}


## Merge with Apt Info----
sub_protein_apt_info<-subset(protein_apt_info,AptName %in% feat_int)
gene_int<-sub_protein_apt_info[match(feat_int,protein_apt_info$AptName),]$EntrezGeneSymbol

predicted_corr_bmi_ptb$Gene_Name = gene_int

## Take Diff Between real vs actual
predicted_corr_bmi_ptb$difference =  est_cov_bmi_ptb - (predicted_corr_bmi_ptb$est_corr_val)

# Plot (sPTB)----

## Cutoffs determined from Permutation Analysis done by Luka 
FDR_cutoff_2 = -0.07538826126066858
FDR_cutoff_1 = 0.07780646670089418
high_diff_proteins<-subset(predicted_corr_bmi_ptb,difference > FDR_cutoff_1)$Gene_Name
others_to_keep = c("IL36A","LAS2","IGSF3","FABP4","EEF2K","CD2","ATP5PB")

all_together = c(high_diff_proteins,others_to_keep)
plot_data <- subset(predicted_corr_bmi_ptb, Gene_Name %in% all_together)

# Reverse sort factor levels alphabetically
plot_data$Gene_Name <- factor(plot_data$Gene_Name, levels = sort(unique(plot_data$Gene_Name), decreasing = TRUE))

plot_data_sPTB = plot_data
bar_plot_pert_sPTB<-ggplot(plot_data,aes(x=Gene_Name,y=difference))+
  geom_col(stat="identity",width = 0.5,color = "black",fill= "gray")+
  coord_flip()+
  #geom_text(aes(label = gene_label),hjust=-0.2,size=2,vjust=0.1)+
  theme_classic()+
  ggtitle("BMI-sPTB diff. corr.\npost-perturb")+
  
  #adsd lines
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = FDR_cutoff_1,linetype = "dashed",color = "red")+
  geom_hline(yintercept = FDR_cutoff_2,linetype = "dashed",color = "red")+
  xlab("Perturbed Protein")+
  ylab("Diff. of Unperturb - Perturb Corr\nBMI to PTB")+
  scale_y_continuous(
    breaks = seq(-0.15,0.2,by=0.1),
    limits = c(-0.2,0.2)
  )+
  theme(#axis.text.y = element_blank(),
    plot.title = element_text(size=16,hjust=0.5,face="bold"),
    axis.ticks.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.title = element_text(face="bold",size=12))

# Save Plots----
ggsave(plot=bar_plot_pert_mPTB,
       path ="../final_plots/figure6",
       filename="Fig6A_left_mPTB_Perturb.pdf",
       units = "in",width=4,height=5.5)

ggsave(plot=bar_plot_pert_sPTB,
       path ="../final_plots/figure6",
       filename="Fig6A_left_sPTB_Perturb.pdf",
       units = "in",width=4,height=5.5)

# Save CSV----
write.csv(plot_data_mPTB[,c(3,4)],"../../Supplemental_data_revisions/fig6a_perturb_mPTB.csv")
write.csv(plot_data_sPTB[,c(3,4)],"../../Supplemental_data_revisions/fig6a_perturb_sPTB.csv")

