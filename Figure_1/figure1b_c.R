# Purpose of Script ---- 
# Figures for Figure 1B and Supplemental Figure 1A

# Load Libraries ----
set.seed(1001)
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggsci)
  library(tidyr)
  library(dplyr)
  library(rstatix)
  library(ggbeeswarm)
  library(ggpubr) #stat_compare_means
  library(smplot2) #for sm_statCorr
})

# linear modeling of PCA----

## Prepare DF ----
patient_info_levels = readRDS("../output_RDS_objs/merged_metadata_and_protein_V1.rds")
proteins_only<-colnames(patient_info_levels)[grepl("seq.",colnames(patient_info_levels))] ## proteins aptamers start with 'seq.'
proteins_only_df<-patient_info_levels[,proteins_only]
proteins_only_df<-as.data.frame(proteins_only_df)
meta_only<-patient_info_levels[,!(colnames(patient_info_levels) %in% proteins_only)]

## Run PCA and merge with meta ----

pca_results = prcomp((proteins_only_df),scale=TRUE)
meta_only<-patient_info_levels[,!(colnames(patient_info_levels) %in% proteins_only)]
merged_meta_pca<-cbind(meta_only,as.data.frame(pca_results$x))

## Linear modeling of PC components and metadata ----
var_interest<-c("bmi30","ga35","sex","study_enrolled","Type","ga_blooddraw","ga_exact","ageatdel","chtn","nullip","ga_blooddraw*Type_mptb","ga_blooddraw*Type_sptb","sex*Type_male","sex*Type_female","age_del*Type_mptb","age_del*Type_sptb")
#var_interest<-c("bmi30","ga35","sex","Type","ga_blooddraw","ga_exact","ageatdel")
PC_num=50
summary_stats_pca<-data.frame(matrix(nrow=length(var_interest)*PC_num,ncol=5))
colnames(summary_stats_pca)<-c("PC_comp","Variable","T_val","p_val","ast")
summary_stats_pca$PC_comp<-paste0("PC",rep(1:PC_num,each=length(var_interest)))
summary_stats_pca$Variable<-rep(var_interest,PC_num)


for (i in unique(summary_stats_pca$PC_comp)){
  sub_PC_df<-summary_stats_pca[summary_stats_pca$PC_comp==i,]
  
  for (j in var_interest){
    sub_eq<- as.formula(paste0("`", i, "` ~","`",j,"`"))
    
    
    if(j=="ga_blooddraw*Type_mptb"){
      sub_sub_pca<-subset(merged_meta_pca,Type %in% c("ctrl","mtpb"))
     sub_eq<- as.formula(paste0("`", i,"`","~","ga_blooddraw","*","Type"))
      lm_model<-lm(sub_eq,data=sub_sub_pca)
      sum_model<-summary(lm_model)
      sum_model_coeff<-sum_model$coefficients
    
      sub_sub_df<-sub_PC_df[sub_PC_df$Variable==j,]
      sub_sub_df$T_val<-sum_model$coefficients[4,3]
      sub_sub_df$p_val<-sum_model$coefficients[4,4]
      
    }
    else if(j=="ga_blooddraw*Type_sptb"){
      sub_sub_pca<-subset(merged_meta_pca,Type %in% c("ctrl","meta_sptb"))
     sub_eq<- as.formula(paste0("`", i,"`","~","ga_blooddraw","*","Type"))
      #sub_eq<- as.formula(paste0("`", i, "` ~","`",j,"`"))
      lm_model<-lm(sub_eq,data=sub_sub_pca)
      sum_model<-summary(lm_model)
      sum_model_coeff<-sum_model$coefficients
    
      sub_sub_df<-sub_PC_df[sub_PC_df$Variable==j,]
      sub_sub_df$T_val<-sum_model$coefficients[4,3]
      sub_sub_df$p_val<-sum_model$coefficients[4,4]
    }
    else if(j=="age_del*Type_mptb"){
      sub_sub_pca<-subset(merged_meta_pca,Type %in% c("ctrl","mtpb"))
     sub_eq<- as.formula(paste0("`", i,"`","~","ageatdel","*","Type"))
      #sub_eq<- as.formula(paste0("`", i, "` ~","`",j,"`"))
      lm_model<-lm(sub_eq,data=sub_sub_pca)
      sum_model<-summary(lm_model)
      sum_model_coeff<-sum_model$coefficients
    
      sub_sub_df<-sub_PC_df[sub_PC_df$Variable==j,]
      sub_sub_df$T_val<-sum_model$coefficients[4,3]
      sub_sub_df$p_val<-sum_model$coefficients[4,4]
    }
    else if(j=="age_del*Type_sptb"){
      sub_sub_pca<-subset(merged_meta_pca,Type %in% c("ctrl","meta_sptb"))
     sub_eq<- as.formula(paste0("`", i,"`","~","ageatdel","*","Type"))
      #sub_eq<- as.formula(paste0("`", i, "` ~","`",j,"`"))
      lm_model<-lm(sub_eq,data=sub_sub_pca)
      sum_model<-summary(lm_model)
      sum_model_coeff<-sum_model$coefficients
    
      sub_sub_df<-sub_PC_df[sub_PC_df$Variable==j,]
      sub_sub_df$T_val<-sum_model$coefficients[4,3]
      sub_sub_df$p_val<-sum_model$coefficients[4,4]
    }
    else if(j=="sex*Type_male"){
       sub_sub_pca<-subset(merged_meta_pca,sex %in% 'Male')
       sub_eq<- as.formula(paste0("`", i,"`","~","Type"))
      lm_model<-kruskal_test(sub_eq,data=sub_sub_pca)
      
      sub_sub_df<-sub_PC_df[sub_PC_df$Variable==j,]
      sub_sub_df$T_val<-lm_model$statistic
      sub_sub_df$p_val<-lm_model$p
    }
    
    else if(j=="sex*Type_female"){
       sub_sub_pca<-subset(merged_meta_pca,sex %in% 'Female')
       sub_eq<- as.formula(paste0("`", i,"`","~","Type"))
      lm_model<-kruskal_test(sub_eq,data=sub_sub_pca)

      sub_sub_df<-sub_PC_df[sub_PC_df$Variable==j,]
      sub_sub_df$T_val<-lm_model$statistic
      sub_sub_df$p_val<-lm_model$p
    }
    else if(class(merged_meta_pca[,j])=="character"){
      if( length(unique(merged_meta_pca[,j])) > 2){
        lm_model<-kruskal_test(sub_eq,data=merged_meta_pca)
      }
      else{
        lm_model<-wilcox_test(sub_eq,data=merged_meta_pca)
      }
      sub_sub_df<-sub_PC_df[sub_PC_df$Variable==j,]
      sub_sub_df$T_val<-lm_model$statistic
      sub_sub_df$p_val<-lm_model$p
      
    }
    
    else{
       lm_model<-lm(sub_eq,data=merged_meta_pca)
       sum_model<-summary(lm_model)
      sum_model_coeff<-sum_model$coefficients
    
      sub_sub_df<-sub_PC_df[sub_PC_df$Variable==j,]
      sub_sub_df$T_val<-sum_model$coefficients[2,3]
      sub_sub_df$p_val<-sum_model$coefficients[2,4]
    }
    
    sub_PC_df[sub_PC_df$Variable==j,]<-sub_sub_df
    
  }
  summary_stats_pca[summary_stats_pca$PC_comp==i,]<-sub_PC_df
}

#addins asterisks 
Signif <- symnum(summary_stats_pca$p_val, corr = FALSE, na = FALSE, cutpoints = c(0, 
    0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))
summary_stats_pca$ast<-Signif

# (Fig 1b) Setup and Plotting PCA Heatmap ----
summary_stats_pca$PC_comp<-factor(summary_stats_pca$PC_comp,levels=unique(summary_stats_pca$PC_comp))
summary_stats_pca$Variable<-factor(summary_stats_pca$Variable,levels=unique(summary_stats_pca$Variable[order(summary_stats_pca$Variable,decreasing = TRUE)]))

## Only plottign first ten and variable of interest ----
sub_PC_int<-paste0("PC",rep(1:10))
sub_summary_stats_PCA<-subset(summary_stats_pca, PC_comp %in% sub_PC_int)
main_var<-c("bmi30","ga35","ga_blooddraw","Type","sex","ageatdel","ga_exact")
sub_summary_stats_PCA_MAIN<-subset(sub_summary_stats_PCA,Variable %in% main_var)
sub_summary_stats_PCA_MAIN$Variable<-factor(sub_summary_stats_PCA_MAIN$Variable,levels=unique(sub_summary_stats_PCA_MAIN$Variable))
sub_summary_stats_PCA_MAIN$Variable<-factor(sub_summary_stats_PCA_MAIN$Variable,labels=c("<30 vs >=30 BMI (kg/m2)","Delivered Before 35 weeks GA?","Fetal Sex","mPTB or sPTB or Ctrl", "GA at Blooddraw","GA of Delivery","Maternal Age"))
sub_summary_stats_PCA_MAIN$PC_comp<-factor(sub_summary_stats_PCA_MAIN$PC_comp,levels=paste0("PC",rep(1:10)))

## Plotting of PCA Heatmap ----
p<-ggplot(sub_summary_stats_PCA_MAIN,aes(x=PC_comp,y=Variable,fill=-log10(p_val)))
heatmap_PCA_MAIN<-p+geom_tile(colour="gray")+
  scale_fill_gradient2(low="#1207A3",mid="#FFFEFE",high="#BB0103",midpoint=-log10(0.05))+
  geom_text(aes(label=ast),color="black",size=6)+
  theme(panel.border=element_rect(colour="black",fill=NA),panel.background = element_blank(),
        axis.title=element_blank(),
        axis.text.x = element_text(size=12,angle=45,hjust=1,color="black"),
        axis.text.y=element_text(size=9,colour="black"),
        legend.title=element_text(size=10),
        legend.text = element_text(size=8))

# (Sup Fig 1A) Plotting of Variance explained
eigs <- pca_results$sdev^2
proportions<-eigs[1:10] / sum(eigs)
PC_list<-paste0("PC",rep(1:10))
df_PC_proportions<-as.data.frame(cbind(PC_list,proportions))
df_PC_proportions$proportions<-as.numeric(df_PC_proportions$proportions)

df_PC_proportions$PC_list<-factor(df_PC_proportions$PC_list,levels=c(PC_list))

plot_pca_var<-ggplot(df_PC_proportions,aes(x=PC_list,y=round(proportions*100,2),group=1))+
  geom_bar(stat="identity")+
  geom_point()+
  geom_line(color="red")+
  ylab("Explained Variation (%)")+
  xlab("Principal component (PC)")+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1))


# (Fig 1C)

## Setup of Colors Used
val_cols<-c("Control"="#868686FF","sPTB"="#0073C2FF","mPTB"="#CD534CFF")
bmi_cols<-pal_npg()(2)
## Setup PC6----
merged_meta_pca$Type = factor(merged_meta_pca$Type,levels=c("ctrl","meta_sptb","mtpb"))
merged_meta_pca$Type = factor(merged_meta_pca$Type,labels=c("Control","sPTB","mPTB"))
PC6_var<-round(df_PC_proportions[df_PC_proportions$PC_list=="PC6",]$proportions*100,2)
my_comparisons <- list( c("Control", "sPTB"), c("sPTB", "mPTB"),c("Control", "mPTB") )


## PCA6 vs BMI ----
bmi_cols<-pal_npg()(2)

PC6_BMI_plot<-ggplot(merged_meta_pca,aes(x=bmi30,y=PC6))+

  geom_violin(aes(fill=bmi30),draw_quantiles=c(0.25,0.5,0.75),alpha=0.3)+
  geom_beeswarm(aes(fill=bmi30),pch=21,size=4,cex=2)+
  scale_fill_manual(guide="none",values=bmi_cols)+
  stat_compare_means(comparisons = list(c("<30", ">=30")), 
                     label = "p.signif",label.x=1.45,size=4,tip.length = 0.01,na.rm=TRUE,bracket.size=0.5,vjust=0.4)+
  xlab("BMI (kg/m2)")+
  ylab(paste0("PC6"," (",PC6_var,"%)"))+
  theme_classic()+
  theme(axis.title=element_text(size=15,face="bold"),
        axis.text=element_text(size=12),
        axis.line=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

## PCA6 vs Type ----

stat.test <- merged_meta_pca %>%
  dunn_test(as.formula(paste0("PC6 ~","Type")), p.adjust.method = "fdr",detailed=TRUE) %>%
  add_xy_position(x = "Type") %>%
  remove_ns()

PC6_Type_plot<-ggplot(merged_meta_pca,aes(x=Type,y=PC6))+
  geom_violin(aes(fill=Type),draw_quantiles=c(0.25,0.5,0.75),alpha=0.3)+
  geom_beeswarm(aes(fill=Type),pch=21,size=4,cex=2)+
  scale_fill_manual(guide="none",name="",values=val_cols)+

  stat_pvalue_manual(stat.test,   label = "p.adj.signif", tip.length = 0.01,size=4,bracket.size=0.5)+

  ylab(paste0("PC6"," (",PC6_var,"%)"))+
  theme_classic()+
  theme(axis.title.y =element_text(size=15,face="bold"),
        strip.text = element_text(size=17,face="bold"),
        strip.background = element_blank(),
        axis.text=element_text(size=12),
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        axis.title.x=element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_text(face="bold"))

# Sup Figure 1B
PC2_var<-round(df_PC_proportions[df_PC_proportions$PC_list=="PC2",]$proportions*100,2)

PC2_age_plot<-ggplot(merged_meta_pca,aes(x=ageatdel,y=PC2))+
  geom_point(aes(fill=Type),pch=21,size=4)+
  scale_fill_manual(guide="none",values=val_cols)+
  ylab(paste0("PC2"," (",PC2_var,"%)"))+
  xlab("Maternal Age at Delivery (years)")+
  sm_statCorr(label_x=35,label_y=-60,legends=TRUE,borders=TRUE,color="black",linetype="dashed",corr_method = "spearman")+
  theme_classic()+
  theme(axis.title=element_text(size=15,face="bold"),
        axis.text=element_text(size=12),
        axis.line=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# Sup Figure 1C
PC5_var<-round(df_PC_proportions[df_PC_proportions$PC_list=="PC5",]$proportions*100,2)
PC5_age_plot<-ggplot(merged_meta_pca,aes(x=ga_blooddraw,y=PC5))+
  geom_point(aes(fill=Type),pch=21,size=4)+
  scale_fill_manual(guide="none",values=val_cols)+
  ylab(paste0("PC5"," (",PC5_var,"%)"))+
  xlab("GA at blood draw (Weeks)")+
  sm_statCorr(label_x=20,label_y=-50,legends=TRUE,borders=TRUE,color="black",linetype="dashed",corr_method = "spearman")+
  theme_classic()+
  theme(axis.title=element_text(size=15,face="bold"),
        axis.text=element_text(size=12),
        axis.line=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
# Saving Data and Plot ----

## Save the actual data used to plot figure 1b----
write.csv(sub_summary_stats_PCA_MAIN,"../Supplemental_data_revisions/fig1b_heatmap.csv")

## Save the data used to plot figure 1c ----
PC6_BMI = data.frame(Patient_ID = paste0("Patient_",1:nrow(merged_meta_pca)),
                     Group = merged_meta_pca$Type,
                     BMI_group = merged_meta_pca$bmi30,
                     PC6_value = merged_meta_pca$PC6
)

write.csv(PC6_BMI,"../Supplemental_data_revisions/fig1c_pc6_valus.csv")

## save the figure 1b ----
ggsave(filename = "fig1b_PCA_heatmap.pdf",heatmap_PCA_MAIN,units="in",
       path="../../final_plots/figure1",
       width=8,
       height=2.5)


## save figure 1c (right)----
ggsave(filename = "fig1c_PC6_BMI.pdf",
       PC6_BMI_plot,
       units="in",
       path="../../final_plots/figure1",
       width=4,
       height=4)

## save figure 1c (left) ----
ggsave(filename = "fig1c_PC6_Type.pdf",
       PC6_Type_plot,
       units="in",
       path="../../final_plots/figure1",
       width=4,
       height=4)



## save the sup figure 1a ----
ggsave(filename = "sup_fig1a_var_explained.pdf",plot_pca_var,units="in",
       path="../final_plots/sup_figure1",
       width=5,
       height=4)

## save the sup figure 1b ----
ggsave(filename = "sup_fig1b_PC2_vs_age.pdf",PC2_age_plot,units="in",
       path="../final_plots/sup_figure1",
       width=5,
       height=4)


## save the sup figure 1c ----
ggsave(filename = "sup_fig1c_PC5_vs_blooddraw.pdf",PC5_age_plot,units="in",
       path="../final_plots/sup_figure1",
       width=5,
       height=4)
