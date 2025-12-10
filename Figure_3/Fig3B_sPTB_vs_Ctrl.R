# Purpose of Script ---- 
# Figures for Figure 2A-B and D

# Preparing Packages and Functions----
set.seed(1010)
source("../helpful_functions/Systems_Serology_PAL.R")

suppressPackageWarnings({
  library(ggpubr)
  library(dplyr)
})

# linear modeling of PCA----

## Prepare DF ----
patient_info_levels = readRDS("../output_RDS_objs/merged_metadata_and_protein_V1.rds")
proteins<-colnames(patient_info_levels)[grepl("seq.",colnames(patient_info_levels))]
merged_data_mod<-patient_info_levels

## Fix some of the variable names ---
merged_data_mod$Type<-gsub("meta_sptb","sptb",merged_data_mod$Type)
merged_data_mod$Type<-gsub("mtpb","mptb",merged_data_mod$Type)

merged_data_mod$bmi30<-gsub("<30","less30",merged_data_mod$bmi30)
merged_data_mod$bmi30<-gsub(">=30","greater30",merged_data_mod$bmi30)
merged_data_mod<-subset(merged_data_mod, Type %in% c("ctrl","sptb"))

## Prepare the X for PLSDA 
rownames(merged_data_mod)<-paste0("studyid_",merged_data_mod$studyid)
X_unfilt_data<-merged_data_mod[,grepl("seq.",colnames(merged_data_mod))]
X_scale<-scale(X_unfilt_data,center=TRUE,scale=TRUE)

## Prepare the Y for PLSDA
y_var<-merged_data_mod$Type
y_var<-factor(y_var,levels=c("ctrl","sptb"))

## Work on Apt DF 
protein_apt_info = readRDS("../output_RDS_objs/Somalogic_protein_apt_info.rds")
max_entrez<-protein_apt_info
max_entrez<-max_entrez %>%
  group_by(Single_Entrez) %>%
  slice_max(f_stat) %>%
  ungroup()
protein_apt_info$max_f_stat_single_entrez="N"
protein_apt_info[protein_apt_info$AptName %in% max_entrez$AptName,]$max_f_stat_single_entrez="Y"

## Prepare df containing info of PLSDA plot ---
my_colors=list(Status=c("ctrl"="#868686FF","sptb"="#0073C2FF"))

sub_gene<-subset(protein_apt_info,AptName %in% colnames(X_scale))
sub_gene<-sub_gene[match(colnames(X_scale),sub_gene$AptName),]
df_features <- data.frame(name = colnames(X_scale))
df_features$label <- c(sub_gene$EntrezGeneSymbol)

opts_plot <- list(df_features = df_features,
                  loading_alpha = 1, # transparency for the loadings
                  score_alpha = 1, # transparency for the scores
                  LV_ind = c(1,2), # which LVs to plot
                  color_features = "antigen", # according to which property (defined in df_features) the features should be color-coded
                  colors = my_colors,
                  y_name = "Status",
                  alpha=1,
                  size=4.5,
                  stroke=0.5) 


# Running Model ----

## Get optimum alpha ----
model <- caret::train(x = data.matrix(X_scale), y = y_var,method="glmnet",
                      trControl = trainControl("LOOCV"),tuneLength = 10)
alpha_tuned <- model$bestTune$alpha

## Lasso Regularization----
opts_sel <- list(n_trials = 1000, threshold = 0.60, return_count = TRUE,alpha=alpha_tuned)
sel_features_df <- select_repeat(X_scale,y_var,
                                 selector = select_lasso,
                                 options = opts_sel)

sel_features_count<-as.data.frame(sel_features_df$feature_count)

#keep those that appear 75% of the time 
nonzero_sel_features<-rownames(sel_features_count)[sel_features_count[,1]>600]
sel_features=nonzero_sel_features


# Get Plots of The PLSDA results----
## Figure 3A (top)----
opts_model<-list(n_LV=2)

X_sel<-X_scale[,sel_features] # working with matrix that has selected features
model_sel<-train_ropls(X_sel,y_var,options=opts_model)
opts_plot$LV_ind<-c(1,2)
plsda_biplot<-visualize_ropls_scores(model_sel,y_var,options=opts_plot)
plsda_biplot<-plsda_biplot+theme(axis.title = element_text(size = 15,face="bold"),
                                 axis.text=element_text(size=12))

## Figure 3A (bottom) ----
plt_VIP_plsda<-VIP_plot(model_sel,X_sel,y_var,my_colors)

## Sup Fig 3A----
opts_plot$X <- X_sel
opts_plot$y <-y_var
opts_plot$LV_ind <- 1
opts_plot$mark_enrichment <- TRUE
plt_loadings_bar <- visualize_ropls_loadings_bar(temp_model_sel, options = opts_plot)
plt_loadings_bar_LV1<-plt_loadings_bar+theme(axis.text.x=element_text(size=8,face="bold"), axis.text.y=element_text(size=6,face="bold"),axis.title.x=element_text(size=10))


# Model Validation ----

## Setup Model parameters----
opts = list(n_folds = 5, pt_trials =100, rf_trials=100,save = TRUE, compare_pred = "y")
vals_trial<-list()


method = list(select=feature_select_iter,
              train = train_ropls,
              predict = predict_ropls,
              score = score_accuracy)
## Run Validation----

vals <- cross_validation_unpaired(X_scale,y_var, method, opts,n_trials=10)


## Plot results (Sup Fig 3C Top) ----
plt_validation<- visualize_validate(vals, options = list(y_label = "accuracy")) +
  theme(axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 12),
        axis.text.x=element_text(size=18))

# Plot LV1 Corr with Proteins----

## Prepare DF----

sub_data_scaled_df<-as.data.frame(X_scale)
sub_data_scaled_df$y_var<-merged_data_mod$Type
sub_data_scaled_df$LV1_score<-model_sel@scoreMN[,1]
sub_data_scaled_df$LV2_score<-model_sel@scoreMN[,2]

##Calculate Correlation----
corr_mat_spear<-as.data.frame(matrix(nrow=length(proteins),ncol=4))
colnames(corr_mat_spear)<-c("protein","pval","padj","rho")
corr_mat_spear$protein<-proteins

for(i in proteins){
  val<-sub_data_scaled_df[,i]
  score<-sub_data_scaled_df$LV1_score
  
  
  temp_df<-as.data.frame(cbind(val,score))
  temp_df$val<-as.numeric(temp_df$val)
  temp_df$score<-as.numeric(temp_df$score)
  temp_rho<-cor.test(x=val,y=score,data=temp_df)
  # 
  corr_mat_spear[corr_mat_spear$protein==i,]$rho= temp_rho$estimate
  corr_mat_spear[corr_mat_spear$protein==i,]$pval=temp_rho$p.value
}

corr_mat_spear$padj<-p.adjust(corr_mat_spear$pval,method="fdr")

corr_mat_spear$genes<-subset(protein_apt_info,AptName %in% corr_mat_spear$protein)$Single_gene
corr_mat_spear$proteins<-subset(protein_apt_info,AptName %in% corr_mat_spear$protein)$Single_UniPro
corr_mat_spear$full<-subset(protein_apt_info, AptName %in% corr_mat_spear$protein)$EntrezGeneSymbol
corr_mat_spear$target<-subset(protein_apt_info,AptName %in% corr_mat_spear$protein)$Target
corr_mat_spear$fscore<-subset(protein_apt_info,AptName %in% corr_mat_spear$protein)$f_stat

corr_mat_spear$entrez<-subset(protein_apt_info,AptName %in% corr_mat_spear$protein)$Single_Entrez

corr_mat_spear$keep_Y_N<-subset(protein_apt_info,AptName %in% corr_mat_spear$protein)$max_f_stat_single_entrez

## Prepare Df for Heatmap---
top_feat_num<-25
corr_mat_spear<-corr_mat_spear[order(corr_mat_spear$pval),]
top_feat<-corr_mat_spear$protein[1:top_feat_num]

new_scaled_df<-sub_data_scaled_df[order(sub_data_scaled_df$LV1_score),]
top_feat_heatmap<-new_scaled_df[,which(colnames(new_scaled_df) %in% top_feat)]
top_feat_heatmap<-t(top_feat_heatmap)

## Setup Annotation info----
y_col<-colorRamp2(c(min(new_scaled_df$LV1_score),max(new_scaled_df$LV1_score)),c("#fefdfa","#b7307a"))
row_ha = HeatmapAnnotation(LV1 = new_scaled_df$LV1_score,
                           Status=new_scaled_df$Type,
                           col = list(Status=c("sptb"="#0073C2FF","ctrl"="#868686FF"),
                                      LV1 = y_col),
                           annotation_legend_param=list(
                             LV1 = list(direction="horizontal",border="black",title_position = "topcenter",at=c(-8,0,8),legend_height=unit(2,"cm"),legend_width=unit(4,"cm")),
                             Status = list(nrow = 1,border="black",title_position = "topcenter",plot=FALSE)),
                           show_legend = c(TRUE,FALSE),
                           border=TRUE,
                           show_annotation_name = FALSE)


z_score_col= colorRamp2(c(min(top_feat_heatmap), 0, max(top_feat_heatmap)), c("#1207A3", "#FFFEFE", "#BB0103"))


#get new row labels
temp_sub<-subset(protein_apt_info, AptName %in% rownames(top_feat_heatmap))
new_row_labels<-temp_sub[match(rownames(top_feat_heatmap), temp_sub$AptName),]$EntrezGeneSymbol

## Plot Heatmap----
temp_ht<-Heatmap((as.matrix(top_feat_heatmap)),col=z_score_col,
                 name = "Z-score",
                 top_annotation = row_ha,
                 border_gp=gpar(col="black",lwd=1),
                 
                 width = ncol(top_feat_heatmap)*unit(1 ,"mm"),
                 height = nrow(top_feat_heatmap)*unit(3,"mm"),
                 
                 show_row_dend = FALSE,
                 na_col = "gray",
                 cluster_columns = FALSE,
                 show_column_names = FALSE,
                 cluster_rows = TRUE,
                 row_names_side = "left",
                 show_row_names = TRUE,
                 heatmap_legend_param = list(border="black",
                                             direction="horizontal",
                                             title_position="topcenter",
                                             legend_height=unit(2,"cm"),legend_width=unit(4,"cm")),
                 row_labels=new_row_labels,
                 
                 row_names_gp = grid::gpar(fontsize=8),
                 show_heatmap_legend = TRUE
)

# Save PDFs----

## Fig3A PLSDA
ggsave(filename = "Fig3B_PLSDA_sPTB_biplot.pdf",plt_scores_sel,units="in",
       path="../final_plots/figure3",
       width=4,
       height=4)


ggsave(filename = "Fig3B_PLSDA_sPTB_VIP_plot.pdf",plt_loadings_bar,units="in",
       path="../final_plots/figure3",
       width=4,
       height=2.5)

## Fig3A Heatmap
pdf(file="../final_plots/figure3/Fig3B_LV1_sPTB_proteins_heatmap",height=5,width=4.8)
draw(temp_ht, heatmap_legend_side = "bottom",
     annotation_legend_side="top")
dev.off()

#supp figures
ggsave(filename = "Sup_Fig3B_sPTB_LV1_loadings.pdf",plt_VIP_temp,units="in",
       path="../final_plots/Sup_figure3",
       width=4,
       height=2.5)

ggsave(filename = "Sup_Fig3C_sPTB_bottom_model_val.pdf",plt_val_p2,units="in",
       path="../final_plots/Sup_figure3",
       width=4,
       height=4)

# Save CSV----
## PLSDA Plot----
fig3b_left = data.frame(Patient_ID = paste0("Patient_",1:nrow(model_sel@scoreMN)),
                        Group = as.vector(y_var),
                        LV1 = as.vector(model_sel@scoreMN[,1]),
                        LV2 = as.vector(model_sel@scoreMN[,2]))


write.csv(fig3b_left,"../../Supplemental_data_revisions/fig3b_plsda_mptb_ctrl.csv")

## Corr Heatmap----
heatmap_save_csv = temp_ht@matrix[row_order(temp_ht),]
rownames(heatmap_save_csv) = new_row_labels[row_order(temp_ht)]
colnames(heatmap_save_csv) = paste0("Patient_",1:ncol(heatmap_save_csv))
df_heatmap = rbind(sub_data_scaled_df$LV1_score,sub_data_scaled_df$Type,heatmap_save_csv)
rownames(df_heatmap)[1:2] = c("LV1","Group")

write.csv(df_heatmap,"../../Supplemental_data_revisions/fig3b_heatmap_right.csv")


# Save RDS----
saveRDS(corr_mat_spear,file="../output_RDS_objs/fig3_PLSDA_sPTB_Ctrl.rds")


