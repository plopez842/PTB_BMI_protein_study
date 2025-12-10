# Purpose of R Script ----
# Make volcano plots for Fig4B - LRT of mPTB vs Ctrl

# Load Packages and Funcs----
suppressPackageStartupMessages({
  library(mgcv)
  library(RhpcBLASctl)
  library(foreach)
  library(doSNOW)
  library(parallelly)
})

source("../helpful_functions/fig4_LRT_parallel_funcs.R")
source("../helpful_functions/fig4_univ_plotting.R")

# Load Variables----

## Patient Data
patient_info_levels = readRDS("../output_RDS_objs/merged_metadata_and_protein_V1.rds")
proteins<-colnames(patient_info_levels)[grepl("seq.",colnames(patient_info_levels))]
merged_data_mod<-patient_info_levels

### Fix some of the variable names ---
merged_data_mod$Type<-gsub("meta_sptb","sptb",merged_data_mod$Type)
merged_data_mod$Type<-gsub("mtpb","mptb",merged_data_mod$Type)

merged_data_mod$bmi30<-gsub("<30","less30",merged_data_mod$bmi30)
merged_data_mod$bmi30<-gsub(">=30","greater30",merged_data_mod$bmi30)
merged_data_mod<-subset(merged_data_mod, Type %in% c("ctrl","mptb"))

sub_ctrl<-subset(merged_data_mod,Type %in% c("ctrl"))
sub_mptb<-subset(merged_data_mod,Type %in% c("mptb"))
sub_merged_data<-rbind(sub_ctrl,sub_mptb)


proteins<-colnames(sub_merged_data)[grepl("seq.",colnames(sub_merged_data))]

## Work on Apt DF 
protein_apt_info = readRDS("../output_RDS_objs/Somalogic_protein_apt_info.rds")
max_entrez<-protein_apt_info
max_entrez<-max_entrez %>%
  group_by(Single_Entrez) %>%
  slice_max(f_stat) %>%
  ungroup()
protein_apt_info$max_f_stat_single_entrez="N"
protein_apt_info[protein_apt_info$AptName %in% max_entrez$AptName,]$max_f_stat_single_entrez="Y"


## Log2FC data----
all_comp_limma<-readRDS(file="../output_RDS_objs/pairwise_comp_volc_plot_list_df.rds")

# LRT Permutation----
#PARALLIZATION - START HERE
cr = availableCores() # number of cores
RhpcBLASctl::blas_set_num_threads(round(cr)) # limit core usage


iterator = proteins # define some iterable object


## par params---------------------------------------------
cl <- makeCluster(cr)
registerDoSNOW(cl)
pb <- txtProgressBar(max = length(iterator), style = 3)
progress <- function(n_) setTxtProgressBar(pb, n_)
opts <- list(progress = progress)


## iterate---------------------------------------------
print('Begin normalization')
LRT_results = foreach(i = 1:length(iterator), .combine=rbind,.packages = c('mgcv'), .export = c('func1'),  
                .verbose = TRUE, .options.snow = opts) %dopar% {
                  x<- proteins[i]# define some variable as afunction of i
                  res<-func1(x,sub_merged_data)
                }

close(pb)
stopCluster(cl)

rownames(LRT_results)<-proteins
colnames(LRT_results)<-"pval"
LRT_results<-as.data.frame(LRT_results)
LRT_results$qval<-p.adjust(LRT_results$pval,method="fdr")

# Work on Plotting Fig 4C----
dis_ctrl_g30<-all_comp_limma$mptbVsctrl_g30
dis_ctrl_g30<-dis_ctrl_g30[order(rownames(dis_ctrl_g30)),]
dis_ctrl_l30<-all_comp_limma$mptbVsctrl_l30
dis_ctrl_l30<-dis_ctrl_l30[order(rownames(dis_ctrl_l30)),]

log2_comp_dis<-as.data.frame(matrix(nrow=length(proteins),ncol=6))
colnames(log2_comp_dis)<-c("protein","log2FC_g30","log2FC_l30","pval_g30","pval_l30","padj_LRT_yes")

log2_comp_dis$protein=rownames(mptb_ctrl_l30)
log2_comp_dis$log2FC_g30=mptb_ctrl_g30$logFC
log2_comp_dis$pval_g30<-mptb_ctrl_g30$P.Value
log2_comp_dis$log2FC_l30=mptb_ctrl_l30$logFC
log2_comp_dis$pval_l30=mptb_ctrl_l30$P.Value
log2_comp_dis$padj_LRT_yes="no"

#diff log2FC
log2_comp_dis$log2_diff<-log2_comp_dis$log2FC_g30-(log2_comp_dis$log2FC_l30)

## Merging LRT Results------
log2_comp_dis=log2_comp_dis[match(rownames(LRT_results),log2_comp_dis$protein),]
log2_comp_dis$pval_LRT_actual<-LRT_results$pval
log2_comp_dis$padj_LRT_actual<-LRT_results$padj

#subset significanat proteins
sig_protein_dis<-subset(LRT_results,qval<0.15)
sig_protein_dis<-rownames(subset(sig_protein_dis,pval<0.005))

log2_comp_dis[log2_comp_dis$protein %in% sig_protein_dis,]$padj_LRT_yes="yes"

#add in label
log2_comp_dis$label=protein_apt_info[match(log2_comp_dis$protein,protein_apt_info$AptName),]$EntrezGeneSymbol
log2_comp_dis[log2_comp_dis$padj_LRT_yes=="no",]$label=NA

#order from leat to greatest
sub_df<-subset(log2_comp_dis,padj_LRT_yes %in% "yes")
keep_labels<-sub_df[!(duplicated(sub_df$label)),]
keep_labels<-keep_labels[order(abs(keep_labels$pval_LRT_actual),decreasing=FALSE),]
keep_top_15<-keep_labels$protein[1:20]
log2_comp_dis[!(log2_comp_dis$protein %in% keep_top_15),]$label=NA

## Plotting volcano plot----

log2_bmi_dis_plot<-ggplot(log2_comp_dis,aes(x=log2FC_l30,y=log2FC_g30,label=label))+
  # Gray points
  geom_point(data = subset(log2_comp_dis, padj_LRT_yes == "no"), 
             aes(color = padj_LRT_yes, alpha = padj_LRT_yes), size = 0.75) +
  # Red points and text
  geom_point(data = subset(log2_comp_dis, padj_LRT_yes == "yes"), 
             aes(color = padj_LRT_yes, alpha = padj_LRT_yes), size = 0.75) +
  geom_text_repel(box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"),
                  force=2,
                  max.overlaps = 15,
                  size=2.5,seed=555)+
  scale_alpha_manual(values=c("yes"=1,"no"=0.5))+
  
  ylim(-1.5,1.5)+
  xlim(-2,2)+
  scale_color_manual(values=c("yes"="#CD534CFF","no"="gray"))+
  
  geom_hline(yintercept=0,linetype="dashed",color="darkgray")+
  geom_vline(xintercept=0,linetype="dashed",color="darkgray")+
  theme_classic()+
  ggtitle("LRT of mPTB vs Ctrl",subtitle=paste0("Total DEP: ",sum(log2_comp_dis$padj_LRT_yes=="yes")))+
  theme(axis.line = element_blank(),
        panel.border=element_rect(colour = "black", fill=NA, size=1),
        axis.title=element_text(face="bold",size=8),
        plot.title=element_text(hjust=0.5,face="bold",size=10),
        plot.subtitle=element_text(hjust=0.5,size=10),
        axis.text = element_text(size=6),
        legend.position = "none")+
  xlab("log2FC(mPTB/Ctrl)\n<30 BMI kg/m3")+
  ylab("log2FC(mPTB/Ctrl)\n>=30 BMI kg/m3")

# Plot Univariates----

##Determine which proteins to plot----
sub_df_feat_int<-log2_comp_dis[!(is.na(log2_comp_dis$label)),]
sub_df_feat_int=(sub_df_feat_int[order(abs(sub_df_feat_int$log2_diff),decreasing=TRUE),])

num_feat = 8
apt_int = sub_df_feat_int[1:num_feat,]
protein_int <- apt_int$protein

## Plot----

color_int = c("#868686FF","#CD534CFF")

stat_test_temp <- list()
for(temp_apt in protein_int){
  stat_test = LRT_univ_fig3_stats(apt_name = temp_apt,
                                  df_int = merged_data_mod)
  stat_test_temp[[temp_apt]] = stat_test
  
}
plot_ind_protein <- list()
for(temp_apt in protein_int){
  stat_test = stat_test_temp[[temp_apt]]
  
  temp_plot = LRT_univ_fig3_plot(temp_apt,
                                 df_int = merged_data_mod,
                                 stat_test = stat_test,
                                 apt_conv = protein_apt_info,
                                 color_groups = color_int)
  plot_ind_protein [[temp_apt]] = temp_plot
  
}

#Save Plot----

## Volcano Plot (4C)
ggsave(log2_bmi_dis_plot,
       path ="../final_plots/figure4",
       filename = "Fig4C_LRT_BMI_mPTB_vs_ctrl.pdf",
       units="in",height=3.5,width=3.5)
## Univariates


lapply(names(plot_ind_protein),
       function(x) ggsave(plot = plot_ind_protein[[x]],
                          filename = paste0("Fig4C_",subset(protein_apt_info, AptName %in% x)$Single_gene,".pdf"),
                          path = "../final_plots/figure4",
                          units = "in",height = 2.5,width=2.5))

# Save CSV files----

##Volcano Plot----
log2_comp_dis$protein_name = protein_apt_info$EntrezGeneSymbol

df_protein_volcano = data.frame(protein = log2_comp_dis$protein_name,
                                xaxis_log2fc_l30 = log2_comp_dis$log2FC_l30,
                                yaxis_log2fc_g30 = log2_comp_dis$log2FC_g30,
                                LRT_padj = log2_comp_dis$pval_LRT_actual,
                                significant = log2_comp_dis$padj_LRT_yes)

write.csv(df_protein_volcano,"../../Supplemental_data_revisions/fig4c_LRT_mPTB_ctrl.csv")

##Univariates----

univeriate_df = data.frame(bmi_group = merged_data_mod$bmi30,
                           PTB_group = merged_data_mod$Type)


protein_df_only =merged_data_mod[,protein_int]
colnames(protein_df_only) = protein_apt_info$EntrezGeneSymbol

univeriate_df = cbind(univeriate_df,protein_df_only)

write.csv(univeriate_df,"../../Supplemental_data_revisions/fig4c_mPTB_univeriates.csv")

# Save RDS (for PCN)----
saveRDS(log2_comp_dis,file="../output_RDS_objs/fig4c_LRT_results_mPTB_Ctrl.rds")

