# Purpose of Script----
# Plot univaraites of proteins that are direct edges to PTB status or BMI
# This is for Figure 5e and Sup Fig 5

# Load Libraries and functions----
source("../helpful_functions/fig5_univ_plotting_func.R")

# Load variables----
df_mPTB = readRDS(file = "../output_RDS_objs/fig5_mPTB_for_PCN.rds")
df_sPTB = readRDS(file = "../output_RDS_objs/fig5_sPTB_for_PCN.rds")

## Work on Apt DF ----
protein_apt_info = readRDS("../output_RDS_objs/Somalogic_protein_apt_info.rds")
max_entrez<-protein_apt_info
max_entrez<-max_entrez %>%
  group_by(Single_Entrez) %>%
  slice_max(f_stat) %>%
  ungroup()
protein_apt_info$max_f_stat_single_entrez="N"
protein_apt_info[protein_apt_info$AptName %in% max_entrez$AptName,]$max_f_stat_single_entrez="Y"


sub_protein_apt_info<-subset(protein_apt_info,AptName %in% PROTEINS_INT)
filtered_protein_apt_info <- protein_apt_info %>%
  group_by(EntrezGeneSymbol) %>%          # Group by the EntrezGeneSymbol column
  slice_max(f_stat, n = 1, with_ties = FALSE) %>%  # Keep the row with the highest f_stat
  ungroup() 

# Work on Main Fig (Figure 5e)----

file_dir = "../final_plots/figure5"


list_of_genes = c("UBE2G2","LRRC25","TMEM132A","LEPR","LAS2","PROK2")

for(gene_to_plot in list_of_genes){
  mPTB_plot<-plotting_univariates(df_mPTB,
                                  PCN_data = mPTB_PCN,
                                  seq_to_gene_df = filtered_protein_apt_info,
                                  gene_int = gene_to_plot ,
                                  PTB_group = "mPTB")+theme(axis.title.x = element_blank())
  
  sPTB_plot<-plotting_univariates(df_sPTB,
                                  PCN_data =sPTB_PCN,
                                  seq_to_gene_df =filtered_protein_apt_info,
                                  gene_int = gene_to_plot ,
                                  PTB_group = "sPTB")+theme(axis.title = element_blank())
  
  plots_merged<-ggarrange(mPTB_plot,sPTB_plot)
  
  plots_merged_annotated <- annotate_figure(
    plots_merged,
    top = text_grob(gene_to_plot , face = "bold", size = 12,hjust=0.4),
    bottom = text_grob("BMI (kg/m3)", size = 10,face="bold",hjust=0.35)
  )
  
  ggsave(plots_merged_annotated,file = paste0(file_dir,gene_to_plot,".pdf"),
         units = "in",width=4,height=2)
}

# Work on Supplemental Figure 5----
list_of_genes = c("CCL22","CELA1","EEF2K","FABP3","LAS2","LCT","LEPR","LRRC25","PROK2","TMEM132A","UBE2G2")
file_dir = "../final_plots/sup_figure5"


for(gene_to_plot in list_of_genes){
  mPTB_plot<-plotting_univariates(df_mPTB,
                                  PCN_data = mPTB_PCN,
                                  seq_to_gene_df = filtered_protein_apt_info,
                                  gene_int = gene_to_plot ,
                                  PTB_group = "mPTB")+theme(axis.title.x = element_blank())
  
  sPTB_plot<-plotting_univariates(df_sPTB,
                                  PCN_data =sPTB_PCN,
                                  seq_to_gene_df =filtered_protein_apt_info,
                                  gene_int = gene_to_plot ,
                                  PTB_group = "sPTB")+theme(axis.title = element_blank())
  
  plots_merged<-ggarrange(mPTB_plot,sPTB_plot)
  
  plots_merged_annotated <- annotate_figure(
    plots_merged,
    top = text_grob(gene_to_plot , face = "bold", size = 12,hjust=0.4),
    bottom = text_grob("BMI (kg/m3)", size = 10,face="bold",hjust=0.35)
  )
  
  ggsave(plots_merged_annotated,file = paste0(file_dir,gene_to_plot,".pdf"),
         units = "in",width=4,height=2)
}



