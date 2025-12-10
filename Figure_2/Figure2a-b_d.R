# Purpose of Script ---- 
# Figures for Figure 2A-B and D

# Preparing Packages and Functions----
set.seed(1010)
source("../helpful_functions/Systems_Serology_PAL.R")
source("../helpful_functions/fig2d_univ_plot.R")

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

#sub_merge_data<-subset(merged_data_mod, Type %in% c("ctrl","mptb", "sptb"))
## Prepare the X for PLSDA 
rownames(merged_data_mod)<-paste0("studyid_",merged_data_mod$studyid)
X_unfilt_data<-merged_data_mod[,grepl("seq.",colnames(merged_data_mod))]
X_scale<-scale(X_unfilt_data,center=TRUE,scale=TRUE)

## Prepare the Y for PLSDA
y_var<-merged_data_mod$Type
y_var<-factor(y_var,levels=c("ctrl","mptb","sptb"))

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
my_colors=list(Status=c("ctrl"="#868686FF","mptb"="#CD534CFF","sptb"="#0073C2FF"))

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
opts_sel <- list(n_trials = 1000, threshold = 0.70, return_count = TRUE,alpha=alpha_tuned)
sel_features_df <- select_repeat(X_scale,y_var,
                                 selector = select_lasso,
                                 options = opts_sel)

sel_features_count<-as.data.frame(sel_features_df$feature_count)

#keep those that appear 70% of the time 
nonzero_sel_features<-rownames(sel_features_count)[sel_features_count[,1]>700]
sel_features=nonzero_sel_features

# Get Plots of The PLSDA results----
## Figure 2A (top)----
opts_model<-list(n_LV=2)

X_sel<-X_scale[,sel_features] # working with matrix that has selected features
model_sel<-train_ropls(X_sel,y_var,options=opts_model)
opts_plot$LV_ind<-c(1,2)
plsda_biplot<-visualize_ropls_scores(model_sel,y_var,options=opts_plot)
plsda_biplot<-plsda_biplot+theme(axis.title = element_text(size = 15,face="bold"),
                                     axis.text=element_text(size=12))

## Figure 2A (bottom) ----
plt_VIP_plsda<-VIP_plot(model_sel,X_sel,y_var,my_colors)

## Figure 2B (LV1)----
opts_plot$X <- X_sel
opts_plot$y <-y_var
opts_plot$LV_ind <- 1
opts_plot$mark_enrichment <- TRUE
plt_loadings_bar <- visualize_ropls_loadings_bar(temp_model_sel, options = opts_plot,order_by = "LV")
plt_loadings_bar_LV1<-plt_loadings_bar+theme(axis.text.x=element_text(size=8,face="bold"), axis.text.y=element_text(size=6,face="bold"),axis.title.x=element_text(size=10))

## Figure 2B (LV2)
opts_plot$X <- X_sel
opts_plot$y <-y_var
opts_plot$LV_ind <- 2
opts_plot$mark_enrichment <- TRUE
plt_loadings_bar <- visualize_ropls_loadings_bar(temp_model_sel, options = opts_plot,order_by = "LV")
plt_loadings_bar_LV2<-plt_loadings_bar+theme(axis.text.x=element_text(size=8,face="bold"), axis.text.y=element_text(size=6,face="bold"),axis.title.x=element_text(size=10))


# Univariate Plots----
## Get top 2 highest and lowest proteins of each loadings

df_loadings <- data.frame(LV1 = ropls::getLoadingMN(model_sel)[,1],
                          LV2 = ropls::getLoadingMN(model_sel)[,2],
                          vip_scores = ropls::getVipVn(model_sel))
df_loadings$features <- rownames(df_loadings)
df_loadings$labels <- options$df_features$label[match(rownames(df_loadings), options$df_features$name)]


# top 2 highest LV1
df_loadings_highest_LV1<-df_loadings[order(df_loadings$LV1,decreasing=TRUE),]
# top 2 lowest LV1
df_loadings_lowest_LV1<-df_loadings[order(df_loadings$LV1,),]

# top 2 highest LV2
df_loadings_highest_LV2<-df_loadings[order(df_loadings$LV2,decreasing=TRUE),]

# top 2 lowest LV2
df_loadings_lowest_LV2<-df_loadings[order(df_loadings$LV2,decreasing),]

#merge the proteins
top_LV_proteins = c(df_loadings_highest_LV1$features[1:2],
                    df_loadings_lowest_LV1$features[1:2],
                    df_loadings_highest_LV2$features[1:2],
                    df_loadings_lowest_LV2$features[1:2])


LV_highest_univariates<-list()
for(protein_int in top_LV_proteins){
  
  LV1_highest_univariates[[protein_int]]<-plotting_3way_PLSDA_univariates(df=merged_data_mod,
                                                                feat=protein_int,
                                                                ref_protein_Apt_mat = protein_apt_info,
                                                                group_colors=my_colors$Status)
}

univ_plots_merged<-ggarrange(plotlist=LV1_highest_univariates,nrow=2,ncol=4)

# Save plots----
## Fig 2A----
gsave(filename = "Fig2A_PLSDA_3way_plot.pdf",plt_scores_sel,units="in",
      path="../final_plots/figure2",
      width=4,
      height=4)

ggsave(filename = "Fig2A_VIP_plot.pdf",plt_VIP_temp,units="in",
       path="../final_plots/figure2",
       width=5,
       height=3)

## Fig 2B----
ggsave(filename = "Fig2B_mPTB_vs_sPTB_ctrl_PLSDA_LV1.pdf",plt_loadings_bar_LV1,units="in",
       path="../final_plots/figure2",
       width=4,
       height=9)

ggsave(filename = "Fig2B_mPTB_vs_sPTB_ctrl_PLSDA_LV2.pdf",plt_loadings_bar_LV2,units="in",
       path="../final_plots/figure2",
       width=4,
       height=9)

## Fig 2D
ggsave(filename = "Fig2D_PLSDA_proteins_univariates.pdf",univ_plots_merged,units="in",
       path="../final_plots/figure2",
       width=10,
       height=5)

# Get the CSV of each plot----
##Fig2A----
df_fig2a = data.frame(Patient_ID = paste0("Patient_",1:nrow(model_sel@scoreMN)),
                      LV1 = model_sel@scoreMN[,1],
                      LV2 = model_sel@scoreMN[,2])

rownames(df_fig2a) = NULL
write.csv(df_fig2a,"../Supplemental_data_revisions/fig2a_REAL_plsda.csv")

##Fig2A (VIP)----
vip_df = model_sel@vipVn[model_sel@vipVn>1]
sub_protein_apt_info = subset(protein_apt_info,AptName %in% names(vip_df))
match(names(vip_df),sub_protein_apt_info$AptName)

vip_df = data.frame(gene_name = sub_protein_apt_info$EntrezGeneSymbol,
                    AptName = names(vip_df),
                    VIP_score = as.vector(vip_df))

vip_df = vip_df[order(vip_df$VIP_score,decreasing=TRUE),]
vip_df$Enriched_in = "mPTB"

vip_df[vip_df$gene_name %in% c("UBASH3A","PHB2","EEF2K","ATP5PB"),]$Enriched_in = "sPTB"

write.csv(vip_df_V2[1:15,c(1,3,4)],"../Supplemental_data_revisions/fig2a_REAL_VIP_plot.csv")


## Fig 2B----

write.csv(df_loadings,"../Supplemental_data_revisions/fig2b_REAL_LV_plots.csv")

##Fig2D----
df_univ = data.frame(Group = merged_data_mod$Type,
                     merged_data_mod[,top_LV_proteins])

rownames(df_univ) = NULL

write.csv(df_univ,"../Supplemental_data_revisions/fig2d_univ.csv")

# Get corr againts LV1----

## Setup Df----
sub_data_scaled_df<-as.data.frame(X_scale)
sub_data_scaled_df$y_var<-merged_data_mod$Type
sub_data_scaled_df$LV1_score<-model_sel@scoreMN[,1]
sub_data_scaled_df$LV2_score<-model_sel@scoreMN[,2]


##Perform Correlation----

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

##Save RDS Objects----
saveRDS(corr_mat_spear,file="../output_RDS_objs/fig2_PLSDA_mPTB_sPTB_Ctrl.rds")


