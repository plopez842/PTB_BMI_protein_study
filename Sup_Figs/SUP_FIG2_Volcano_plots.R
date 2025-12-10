# Purpose of R Script ----
# Make volcano plots for Sup Fig2

# load libraries----

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  
})
# Load variables----

## Prepare DF ----
patient_info_levels = readRDS("../output_RDS_objs/merged_metadata_and_protein_V1.rds")
proteins<-colnames(patient_info_levels)[grepl("seq.",colnames(patient_info_levels))]
merged_data_mod<-patient_info_levels

## Fix some of the variable names ---
merged_data_mod$Type<-gsub("meta_sptb","sptb",merged_data_mod$Type)
merged_data_mod$Type<-gsub("mtpb","mptb",merged_data_mod$Type)

merged_data_mod$bmi30<-gsub("<30","less30",merged_data_mod$bmi30)
merged_data_mod$bmi30<-gsub(">=30","greater30",merged_data_mod$bmi30)
merged_data_mod$group<-paste0(merged_data_mod$Type,"_",merged_data_mod$bmi30)
group<-merged_data_mod$group

## Work on Apt DF 
protein_apt_info = readRDS("../output_RDS_objs/Somalogic_protein_apt_info.rds")
max_entrez<-protein_apt_info
max_entrez<-max_entrez %>%
  group_by(Single_Entrez) %>%
  slice_max(f_stat) %>%
  ungroup()
protein_apt_info$max_f_stat_single_entrez="N"
protein_apt_info[protein_apt_info$AptName %in% max_entrez$AptName,]$max_f_stat_single_entrez="Y"

# Run LM----

## Prepare Contrasts----

data_only<-merged_data_mod[,proteins]
design = model.matrix(~0 + group)

vfit = lmFit(t(data_only), design)


contrasts<-makeContrasts(
  mptbVsctrl=(groupmptb_greater30 + groupmptb_less30)/2 - (groupctrl_less30 + groupctrl_greater30)/2,
  sptbVsctrl=(groupsptb_greater30 + groupsptb_less30)/2 - (groupctrl_less30 + groupctrl_greater30)/2,
  sptbVsmptb=(groupsptb_greater30 + groupsptb_less30)/2 - (groupmptb_greater30 + groupmptb_less30)/2,
  g30Vsl30=(groupctrl_greater30+groupmptb_greater30+groupsptb_greater30)/3 - (groupctrl_less30+groupmptb_less30+groupsptb_less30)/3,
  mptbVsctrl_g30=groupmptb_greater30-groupctrl_greater30,
  sptbVsctrl_g30=groupsptb_greater30-groupctrl_greater30,
  sptbVsmptb_g30=groupsptb_greater30-groupmptb_greater30,
  mptbVsctrl_l30=groupmptb_less30-groupctrl_less30,
  sptbVsctrl_l30=groupsptb_less30-groupctrl_less30,
  sptbVsmptb_l30=groupsptb_less30-groupmptb_less30,
  g30Vsl30_ctrl=groupctrl_greater30-groupctrl_less30,
  g30Vsl30_mptb=groupmptb_greater30-groupmptb_less30,
  g30Vsl30_sptb=groupsptb_greater30-groupsptb_less30,
  levels=colnames(design)
)

fit2 <- contrasts.fit(vfit, contrasts)
efit = eBayes(fit2)

comp<-c("mptbVsctrl","sptbVsctrl","sptbVsmptb","g30Vsl30","mptbVsctrl_g30","sptbVsctrl_g30","sptbVsmptb_g30","mptbVsctrl_l30","sptbVsctrl_l30","sptbVsmptb_l30","g30Vsl30_ctrl","g30Vsl30_mptb","g30Vsl30_sptb")

comp_list<-list()

for(i in 1:length(comp)){
  comp_int<-topTable(efit,coef=i,number=length(proteins))
  
  #add in single gene and full gene
  order_protein_info<-protein_apt_info[match(rownames(comp_int),protein_apt_info$AptName),]
  
  comp_int$Single_gene<-order_protein_info$Single_gene
  comp_int$Single_prot<-order_protein_info$Single_UniPro
  comp_int$EntrezGeneID<-order_protein_info$EntrezGeneSymbol
  
  name_comp<-comp[i]
  
  comp_list[[name_comp]]<-comp_int
  
}

# Volcano Plots----

## Prep DF----
merged_comp<-comp_list

adj_pval_thresh<-0.05
for (i in 1:length(merged_comp)){
  merged_comp[[i]]<-merged_comp[[i]] %>%
    mutate(phenotype = case_when (
      adj.P.Val < adj_pval_thresh~ "yes",
      adj.P.Val > adj_pval_thresh~ "no",
      adj.P.Val < adj_pval_thresh~ "no",
      adj.P.Val > adj_pval_thresh~ "no"
      
    ))
  
}

for (i in 1:length(merged_comp)){
  merged_comp[[i]]$label=NA
  merged_comp[[i]][merged_comp[[i]]$phenotype=="yes",]$label=merged_comp[[i]][merged_comp[[i]]$phenotype=="yes",]$EntrezGeneID
  
  
}
## Plot----
names_volc_plot<-c("mPTB vs Ctrl",
                   "sPTB vs Ctrl",
                   "sPTB vs mPTB",
                   ">=30 vs <30 BMI",
                   
                   ">=30 BMI: mPTB vs Ctrl",
                   ">=30 BMI: sPTB vs Ctrl",
                   ">=30 BMI: sPTB vs mPTB",
                   "<30 BMI: mPTB vs Ctrl",
                   "<30 BMI: sPTB vs Ctrl",
                   "<30 BMI: sPTB vs mPTB",
                   "Ctrl: >=30 vs <30 BMI",
                   "mPTB: >=30 vs <30 BMI",
                   "sPTB: >=30 vs <30 BMI")

volc_plot<-list()
for (i in 1:length(merged_comp)){
  title_name<-names_volc_plot[i]
  temp_df<-merged_comp[[i]]
  temp_plot<-ggplot(temp_df,aes(x=logFC,y=-log10(adj.P.Val),label=label))+
    geom_point(aes(fill=phenotype),pch=21)+
    scale_fill_manual(values=c("no"="black","yes"="red"))+
    theme_classic()+
    ylab("-log10(FDR)")+
    xlab("log2FC")+
    ggtitle(title_name)+
    geom_vline(xintercept=0,linetype="dashed",color="gray")+
    geom_hline(yintercept=1.30103,color="gray")+
    theme(legend.position = "none",
          axis.title=element_text(face="bold",size=8),
          plot.title=element_text(hjust=0.5,face="bold",size=10),
          plot.subtitle = element_text(hjust=0.5,size=9.5),
          axis.text = element_text(size=6))
  volc_plot[[i]]<-temp_plot
}
names(volc_plot)<-names_volc_plot

# Save Plots

ggsave(volc_plot$`mPTB vs Ctrl`,
       path="../../final_plots/Sup_figure2",
       file="Sup_Fig2a_mPTB_Ctrl.pdf",
       units="in",width=2.75,height=2.75)

ggsave(volc_plot$`sPTB vs Ctrl`,
       path="../../final_plots/Sup_figure2",
       file="Sup_Fig2b_sPTB_Ctrl.pdf",
       units="in",width=2.75,height=2.75)

ggsave(volc_plot$`sPTB vs mPTB`,
       path="../../final_plots/Sup_figure2",
       file="Sup_Fig2c_sPTB_mPTB.pdf",
       units="in",width=2.75,height=2.75)

# Save RDS----
saveRDS(merged_comp,file="../output_RDS_objs/pairwise_comp_volc_plot_list_df.rds")

