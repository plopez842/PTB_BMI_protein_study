
# Purpose of Script----
# Get the full and partial matrix from graph lasso etc 
# We also get supplemental figures 4a-d (only for sPTB)

# Load Functions and Libraries----

functions_path_dir="../helpful_functions/fig4_PCN"
r_files <- list.files(path = functions_path_dir, pattern = "\\.R$", full.names = TRUE)
sapply(r_files,source)


## load functions and libraries
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(colorRamp2)
})

# Datas
## Object and Protein
df_int <- readRDS("../output_RDS_objs/fig5_sPTB_for_PCN.rds")
test_rows<-colnames(df_int)
data_matrix=df_int
data_matrix = preprocess_data(data_matrix)
samplecov<-cov(data_matrix[, test_rows], use = "complete.obs") # not scaled 

## Work on Apt DF 
protein_apt_info = readRDS("../output_RDS_objs/Somalogic_protein_apt_info.rds")
max_entrez<-protein_apt_info
max_entrez<-max_entrez %>%
  group_by(Single_Entrez) %>%
  slice_max(f_stat) %>%
  ungroup()
protein_apt_info$max_f_stat_single_entrez="N"
protein_apt_info[protein_apt_info$AptName %in% max_entrez$AptName,]$max_f_stat_single_entrez="Y"

## Load Glasso 
daterun = "20250304"
min = 0.02
max=0.2
interval = 0.01


#load the stability matrix
filename_stab_mat = paste0(daterun,"_sPTB_glasso_stab_mat_lambda_min_",min,"_lambda_max_",max,".Rds")
full_stabmat = readRDS(paste0("../output_RDS_objs/Fig5_PCN_int_objects/",filename_stab_mat))

#load the user properties
filename_user_param = paste0(daterun,"_sPTB_USER_PARAM_lambda_min_",min,"_lambda_max_",max,".RData")
load(paste0("../output_RDS_objs/Fig5_PCN_int_objects/",filename_user_param))
interval = user_param$lambda_interval
colnames(full_stabmat) = seq(min,max,by=interval)

samplecov<-cov(data_matrix[, test_rows], use = "complete.obs") # not scaled 

# use the results from comp_mat to identify best 'cutoff'----
cutoff = 0.90
min=user_param$lambda_min
max=user_param$lambda_max

stabmat=full_stabmat[,as.character(seq(min,max,by=user_param$lambda_interval))]
p<-length(test_rows)
adjmat<-matrix(as.numeric(apply(stabmat, 1, max) > cutoff), p, p) 


# R2 calcs of resutls (Sup Fig4d, sPTB)----

## get estimated covariance matrix
model_prec_mat<-stable_partcorr(adjmat, samplecov)
est_covmat<-model_prec_mat$covarmat

## get values of indirect correlates
adj_mat_zeros = get_off_diagonal(adjmat)
adj_mat_zeros = adj_mat_zeros==0

##real
samplecov_sub<-sub_cov_adj_mat(adjmat,samplecov)
diag_sample <- get_off_diagonal(samplecov_sub)
diag_sample <- diag_sample[adj_mat_zeros]

##est
sub_est_cov_real<-sub_cov_adj_mat(adjmat,est_covmat)
diag_est <- get_off_diagonal(sub_est_cov_real)
diag_est <- diag_est[adj_mat_zeros]


actual_off_diag <- diag_sample
estimated_off_diag <- diag_est

comp_r2_real_act=data.frame(real=actual_off_diag,
                            est=estimated_off_diag)

r2_label <- paste0("R² = ", format(comp_mat[comp_mat$cutoff==0.9 & comp_mat$var=="act",]$R2, digits = 3))
## Plot----
r2_plots= ggplot(comp_r2_real_act,aes(x=real,y=est))+geom_point()+theme_classic()+
  ylab("Predicted")+
  xlab("Observed")+
  geom_hex(bins = 50) +  # Hexbin density plot
  scale_fill_viridis_c() +  # Use a color gradient for density
  ggtitle("Predicted vs Actual\nCorrelations (sPTB vs Ctrl)")+
  theme(plot.title=element_text(hjust=0.5,face="bold",size=12),
        axis.title = element_text(size=9),
        axis.text = element_text(size=7),
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_blank(),  # Remove legend text
        axis.line = element_blank(),  # Remove both x-axis and y-axis lines
        panel.border = element_rect(color = "black", fill = NA, size = 1.5))+
  geom_abline(slope = 1, intercept = 0, color = "red", size = 1)+  # Add red y=x line
  annotate("text", x = min(comp_r2_real_act$real) + 0.1, y = max(comp_r2_real_act$est) - 0.1, 
           label = r2_label, size = 3, hjust = 0, fontface = "bold")


## save plot----
ggsave(r2_plots,
       path ="../final_plots/sup_figure4", 
       filename = "Sup_Fig4D_R2_estimated_indirect_correlates_sPTB.pdf",
       width =3.5 , height = 3.5,
       units = "in")



# Work on the Actual vs Predicted Full/Partial Heatmap----

## Get the Predicted PCN

##now get the partial correlates----
partial_corr_matrix=-model_prec_mat$precmat/sqrt(outer(diag(model_prec_mat$precmat), diag(model_prec_mat$precmat)))
partial_corr_matrix<-get_PC_glasso_matrix(adjmat,partial_corr_matrix)
## convert partial correlates of diagonal elements to 1 FOR AESTHETIC PURPOSES
diag(partial_corr_matrix) = 1

##HEATMAPs----
col_fun = colorRamp2(c(-1, 0, 1), c("#1207A3", "#FFFEFE", "#BB0103"))
real_cov_heatmap<-Heatmap(samplecov,
                          name="rho",
                          column_title = "Empirical Correlation Matrix\nsPTB vs Ctrl",
                          column_title_gp = gpar(fontface="bold"),
                          border_gp = gpar(col = "black"),
                          col=col_fun,
                          show_row_names=FALSE,
                          show_column_names = FALSE,
                          #column_labels = NA,
                          show_row_dend = FALSE,
                          show_column_dend = FALSE)

real_cov_heat_draw = draw(real_cov_heatmap)


# get the partial correlate matrices
partial<-Heatmap(partial_corr_matrix[row_order(real_cov_heat_draw),column_order(real_cov_heat_draw)],
                 border_gp = gpar(col = "black"),
                 name="rho",
                 column_title = "Partial Correlations\nsPTB vs Ctrl",
                 column_title_gp = gpar(fontface="bold"),
                 col=col_fun,
                 cluster_rows = FALSE,
                 cluster_columns = FALSE,
                 show_row_names=FALSE,
                 show_column_names = FALSE,
                 #column_labels = NA,
                 show_row_dend = FALSE,
                 show_column_dend = FALSE)
partial_draw <- draw(partial)


# get the esimated covariance matrix 
estcov_mat_ht<-Heatmap(est_covmat[row_order(real_cov_heat_draw),column_order(real_cov_heat_draw)],
                       name="rho",
                       border_gp = gpar(col = "black"),
                       column_title = "Est. Correlation Matrix\nfrom PCN model (sPTB vs Ctrl)",
                       column_title_gp = gpar(fontface="bold"),
                       col=col_fun,
                       show_row_names=FALSE,
                       show_column_names = FALSE,
                       #column_labels = NA,
                       show_row_dend = FALSE,
                       show_column_dend = FALSE)
# SAVE HEATMAPS----
# save the real 
path ="../final_plots/sup_figure4"
pdf(file = "../final_plots/sup_figure4/SUP_FIG4C_sPTB_REAL_COV.pdf",width =3.5 , height = 3)
draw(real_cov_heatmap)
dev.off()

#save the PCN
pdf(file = "../final_plots/sup_figure4/SUP_FIG4C_sPTB_PCN.pdf",width =3.5 , height = 3)
draw(partial)
dev.off()

#save the estimated correlation matrix
pdf(file = "../final_plots/sup_figure4/SUP_FIG4C_sPTB_EST.pdf",width =3.5 , height = 3)
draw(estcov_mat_ht)
dev.off()

# Save RDS 

feat_int<-test_rows[3:length(test_rows)]
sub_protein_apt<-subset(protein_apt_info,AptName %in% feat_int)
gene_int<-sub_protein_apt[match(feat_int,sub_protein_apt$AptName),]$EntrezGeneSymbol

## replacing the new names 
partial_corr_GENES<-partial_corr_matrix
diag(partial_corr_GENES)<-0
colnames(partial_corr_GENES)<-c("BMI","PTB_status",gene_int)
rownames(partial_corr_GENES)<-c("BMI","PTB_status",gene_int)


full_corr_GENES<-samplecov
diag(full_corr_GENES)<-0
colnames(full_corr_GENES)<-c("BMI","PTB_status",gene_int)
rownames(full_corr_GENES)<-c("BMI","PTB_status",gene_int)


saveRDS(full_corr_GENES,"../output_RDS_objs/Fig5_PCN_int_objects/Fig5_sPTB_sampleCOV.RDS")
saveRDS(partial_corr_GENES,"../output_RDS_objs/Fig5_PCN_int_objects/Fig5_sPTB_PCN.RDS")

# Precion Matrix (for Perturbation)
saveRDS(model_prec_mat,"../output_RDS_objs/Fig5_PCN_int_objects/Fig5_sPTB_precision_matrix.RDS")
