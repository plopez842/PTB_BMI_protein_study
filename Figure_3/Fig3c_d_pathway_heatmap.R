# Purpose of Script ---- 
# Figures for Figure 3C-D

# Load Packages----
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

source("../helpful_functions/fig3_pathway_plot_functions.R")

# Load Variables----
corr_mat_mPTB = readRDS("../output_RDS_objs/fig3_PLSDA_mPTB_Ctrl.rds")
corr_mat_sPTB = readRDS("../output_RDS_objs/fig3_PLSDA_sPTB_Ctrl.rds")



## Work on Apt DF 
protein_apt_info = readRDS("../output_RDS_objs/Somalogic_protein_apt_info.rds")
max_entrez<-protein_apt_info
max_entrez<-max_entrez %>%
  group_by(Single_Entrez) %>%
  slice_max(f_stat) %>%
  ungroup()
protein_apt_info$max_f_stat_single_entrez="N"
protein_apt_info[protein_apt_info$AptName %in% max_entrez$AptName,]$max_f_stat_single_entrez="Y"

# Run Pathway Analysis----
rho_thresh<-0.25
## Start with mPTB ----

# keep only a single-measurement of each protein, prioritizing the one with higher F_stat
corr_mat_mPTB<-subset(corr_mat_mPTB,keep_Y_N %in% "Y" )
corr_mat_mPTB<-subset(corr_mat_mPTB,genes %in% names( which(table(corr_mat_mPTB$genes)>1)))

#focus on positive 
corr_mat_mPTB_df<-subset(corr_mat_mPTB,padj<0.25)
corr_mat_mPTB_pos<-subset(corr_mat_mPTB_df,rho>rho_thresh)
corr_proteins_disease_mPTB_pos<-subset(protein_apt_info,AptName %in% corr_mat_mPTB_pos$protein)$Single_Entrez

mPTB_pos_path<- enrichPathway(gene=corr_proteins_disease_mPTB_pos,
                                  pvalueCutoff = 0.15,
                                  qvalueCutoff = 0.3,
                                  organism = "human",
                                  universe=unique(protein_apt_info$Single_Entrez),
                                  minGSSize=10,
                                  maxGSSize=6000,
                                  readable = TRUE
)


#focus on negative
corr_mat_mPTB_neg<-subset(corr_mat_mPTB_df,rho<-rho_thresh)
corr_proteins_disease_mPTB_neg<-subset(protein_apt_info,AptName %in% corr_mat_mPTB_neg$protein)$Single_Entrez

mPTB_neg_path<- enrichPathway(gene=corr_proteins_disease_mPTB_neg,
                              pvalueCutoff = 0.15,
                              qvalueCutoff = 0.3,
                              organism = "human",
                              universe=unique(protein_apt_info$Single_Entrez),
                              minGSSize=10,
                              maxGSSize=6000,
                              readable = TRUE
)

## Work wwith sPTB ----

# keep only a single-measurement of each protein, prioritizing the one with higher F_stat
corr_mat_sPTB<-subset(corr_mat_sPTB,keep_Y_N %in% "Y" )
corr_mat_sPTB<-subset(corr_mat_sPTB,genes %in% names( which(table(corr_mat_sPTB$genes)>1)))

#focus on positive 
corr_mat_sPTB_df<-subset(corr_mat_sPTB,padj<0.25)
corr_mat_sPTB_pos<-subset(corr_mat_sPTB_df,rho>rho_thresh)
corr_proteins_disease_sPTB_pos<-subset(protein_apt_info,AptName %in% corr_mat_sPTB_pos$protein)$Single_Entrez

sPTB_pos_path<- enrichPathway(gene=corr_proteins_disease_sPTB_pos,
                              pvalueCutoff = 0.15,
                              qvalueCutoff = 0.3,
                              organism = "human",
                              universe=unique(protein_apt_info$Single_Entrez),
                              minGSSize=10,
                              maxGSSize=6000,
                              readable = TRUE
)


#focus on negative
corr_mat_sPTB_neg<-subset(corr_mat_sPTB_df,rho<-rho_thresh)
corr_proteins_disease_sPTB_neg<-subset(protein_apt_info,AptName %in% corr_mat_sPTB_neg$protein)$Single_Entrez

sPTB_neg_path<- enrichPathway(gene=corr_proteins_disease_sPTB_neg,
                              pvalueCutoff = 0.15,
                              qvalueCutoff = 0.3,
                              organism = "human",
                              universe=unique(protein_apt_info$Single_Entrez),
                              minGSSize=10,
                              maxGSSize=6000,
                              readable = TRUE
)


# Merge Pathway Results ----

mPTB_vs_ctrl_pos<-mPTB_pos_path@result
mPTB_vs_ctrl_neg<-mPTB_neg_path@result
sPTB_vs_ctrl_neg<-sPTB_neg_path@result
# no results for sPTB_pos_path

pathways_list<-list(mPTB_vs_ctrl_pos,
                    mPTB_vs_ctrl_neg,
                    sPTB_vs_ctrl_neg)

names(pathways_list)<-c("mPTB_pos_path","mPTB_neg_path","sPTB_neg_path")

qvalue_thresh<-0.25
sub_pathways_list<-pathways_list

for(i in names(pathways_list)){
  sub_df<-pathways_list[[i]]
  sub_df<-subset(sub_df,qvalue < qvalue_thresh)
  sub_pathways_list[[i]]<-sub_df
}

# Prep for Plotting ----
keep_pos_mPTB_pathways<-c("RAS processing",
                          "Class I MHC mediated antigen processing & presentation",
                          "Antigen processing: Ubiquitination & Proteasome degradation",
                          "MAP2K and MAPK activation",
                          "Transcriptional regulation by RUNX3",
                          "MAPK family signaling cascades")

keep_neg_mPTB_pathways<-c("Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)",
                          "FGFR2 ligand binding and activation",
                          "TNFs bind their physiological receptors",
                          "Signaling by FGFR in disease",
                          "Cytokine Signaling in Immune system",
                          "Interleukin-6 family signaling")

keep_neg_sPTB_pathways<-c("Signaling by Interleukins",
                          "Cytokine Signaling in Immune system",
                          "Extracellular matrix organization",
                          "Interleukin-4 and Interleukin-13 signaling"
)

keep_pathways<-unique(c(keep_pos_mPTB_pathways,keep_neg_mPTB_pathways,keep_neg_sPTB_pathways))
PTB_types<-c("mPTB","mPTB","sPTB")
directionality_group<-c("Enriched\nrelative to Ctrl",
                        "Suppressed\nrelative to Ctrl",
                        "Suppressed\nrelative to Ctrl")
position_to_columns = 1
pathways_to_plot<-sub_pathways_list


for(i in names(pathways_to_plot)){
  sub_df<-pathways_to_plot[[i]]
  sub_df<-subset(sub_df, Description %in% keep_pathways)
  
  #adding types and directionality
  sub_df$Type <- PTB_types[position_to_columns]
  sub_df$Dir <- directionality_group[position_to_columns]
  
  pathways_to_plot[[i]]<-sub_df
  
  position_to_columns = position_to_columns + 1
}

ALL_pathway_groups_merged<-bind_rows(pathways_to_plot)

#modify pathways name
ALL_pathway_groups_merged[ALL_pathway_groups_merged$Description=="Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)",]$Description = "Regulation of IGF transport and\nuptake by IGFBPs"

ALL_pathway_groups_merged[ALL_pathway_groups_merged$Description=="Antigen processing: Ubiquitination & Proteasome degradation",]$Description = "Antigen processing: Ubiquitination\n& Proteasome degradation"

ALL_pathway_groups_merged[ALL_pathway_groups_merged$Description=="Class I MHC mediated antigen processing & presentation",]$Description = "Class I MHC mediated antigen\nprocessing & presentation"



type_col<-c("mPTB"="#CD534CFF","sPTB"="#0073C2FF")

# Fig 3C Plot ----
sub_mPTB_sPTB_OE_neg_plot<-ggplot(ALL_pathway_groups_merged,aes(x=reorder(Description,Count),y=Count,label=Count,fill=Type))+
  geom_col(aes(fill=Type),position = position_dodge2(width = 0.9, preserve = "single"),color="black") +
  geom_text(position = position_dodge2(width = 0.9, preserve = "single"), angle = 0, vjust=0.2,hjust=-0.2)+
  scale_fill_manual(values=type_col)+
  coord_flip()+
  facet_grid(Dir~.,scales="free",space="free")+
  scale_y_continuous(limits = c(0,100),
                     breaks=c(0,50,100),
                     expand=c(0.0001, 0.05)) +
  theme_classic()+
  theme(axis.title.y=element_blank(),
        axis.line=element_blank(),
        strip.placement = "outside", 
        #axis.text.y=element_text(size=8),
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        axis.text.y=element_text(size=11,face="bold"),
        axis.text.x=element_text(size=8,face="bold"),
        plot.title=element_text(hjust=0.5,face="bold",size=12),
        axis.title.x=element_text(size=10,face="bold"),
        strip.background.y=element_rect(color = NA,  fill=NA), 
        #plot.background = element_rect(color="black",fill=NA),
        #plot.
        legend.title=element_text(size=10,face="bold"),
        legend.text=element_text(size=8),
        strip.text.x = element_text(size=9,face="bold"),
        strip.text.y = element_text(angle=0,size=11,face="bold"),
        strip.background = element_rect(color=NULL),legend.position = "bottom",
        legend.background = element_rect(
          size=0.5, linetype="solid", 
          colour ="black"))+
  guides(size=guide_legend(title.position="top",title.hjust=0.5))


# Start on looking at the log2FC of the genes in each pathway ----

## load the log2FC data----
all_comp_limma<-readRDS(file="../output_RDS_objs/pairwise_comp_volc_plot_list_df.rds")

log2FC_mPTB_vs_Ctrl<-all_comp_limma$mptbVsctrl
log2FC_sPTB_vs_Ctrl<-all_comp_limma$sptbVsctrl

#add in apt
log2FC_mPTB_vs_Ctrl$AptName<-rownames(log2FC_mPTB_vs_Ctrl)
log2FC_sPTB_vs_Ctrl$AptName<-rownames(log2FC_sPTB_vs_Ctrl)

## merge the pathway into list----
merge_pathway_output<-list()

merge_pathway_output[["mPTB_neg"]]<-mPTB_neg_path
merge_pathway_output[["sPTB_neg"]]<-sPTB_neg_path

merge_pathway_output[["mPTB_pos"]]<-mPTB_pos_path

## get the correlation into df----

corr_vals_list<-list()

corr_vals_list[["sPTB_neg"]]<-subset(corr_mat_sPTB,rho<=(-0.25) & padj<0.25)
corr_vals_list[["sPTB_neg"]]<-merge(corr_vals_list[["sPTB_neg"]],log2FC_sPTB_vs_Ctrl,by.x="protein",by.y="AptName")

corr_vals_list[["mPTB_neg"]]<-subset(corr_mat_mPTB,rho<=(-0.25)& padj<0.25)
corr_vals_list[["mPTB_neg"]]<-merge(corr_vals_list[["mPTB_neg"]],log2FC_mPTB_vs_Ctrl,by.x="protein",by.y="AptName")

corr_vals_list[["mPTB_pos"]]<-subset(corr_mat_mPTB,rho>0.25 & padj<0.25)
corr_vals_list[["mPTB_pos"]]<-merge(corr_vals_list[["mPTB_pos"]],log2FC_mPTB_vs_Ctrl,by.x="protein",by.y="AptName")

# Plot the Pathway Heatmaps----
## Heatmap 1----
path_int<-("Class I MHC mediated antigen processing & presentation")
which_group<-"mPTB_pos"

apts_of_interest<-apt_from_pathways(path_int,merge_pathway_output,corr_vals_list,which_group)

#get the logFC from each of the groups
sub_log_mPTB<-subset(log2FC_mPTB_vs_Ctrl, AptName %in% apts_of_interest)
sub_log_mPTB<-sub_log_mPTB[order(sub_log_mPTB$AptName),]
sub_log_sPTB<-subset(log2FC_sPTB_vs_Ctrl, AptName %in% apts_of_interest)
sub_log_sPTB<-sub_log_sPTB[order(sub_log_sPTB$AptName),]

df_heatmap<-data.frame(mPTB=sub_log_mPTB$logFC,
                       sPTB=sub_log_sPTB$logFC)
rownames(df_heatmap)=sub_log_mPTB$EntrezGeneID
heatmap_1<-heatmap_plot(df_heatmap)

## Heatmap 2----
path_int<-("Interleukin-4 and Interleukin-13 signaling")
which_group<-"mPTB_neg"

apts_of_interest_MPTB<-apt_from_pathways(path_int,merge_pathway_output,corr_vals_list,which_group)

which_group<-"sPTB_neg"
apts_of_interest_SPTB<-apt_from_pathways(path_int,merge_pathway_output,corr_vals_list,which_group)

#get the logFC from each of the groups
apts_of_interest<- unique(apts_of_interest_SPTB,apts_of_interest_MPTB)
sub_log_mPTB<-subset(log2FC_mPTB_vs_Ctrl, AptName %in% apts_of_interest)
sub_log_mPTB<-sub_log_mPTB[order(sub_log_mPTB$AptName),]
sub_log_sPTB<-subset(log2FC_sPTB_vs_Ctrl, AptName %in% apts_of_interest)
sub_log_sPTB<-sub_log_sPTB[order(sub_log_sPTB$AptName),]


df_heatmap<-data.frame(mPTB=sub_log_mPTB$logFC,
                       sPTB=sub_log_sPTB$logFC)

rownames(df_heatmap)=sub_log_mPTB$EntrezGeneID


heatmap_2<-heatmap_plot(df_heatmap)

## Heatmap 3----
path_int<-("Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)")
which_group<-"mPTB_neg"

apts_of_interest_MPTB<-apt_from_pathways(path_int,merge_pathway_output,corr_vals_list,which_group)

apts_of_interest<- unique(apts_of_interest_MPTB)
sub_log_mPTB<-subset(log2FC_mPTB_vs_Ctrl, AptName %in% apts_of_interest)
sub_log_mPTB<-sub_log_mPTB[order(sub_log_mPTB$AptName),]
sub_log_sPTB<-subset(log2FC_sPTB_vs_Ctrl, AptName %in% apts_of_interest)
sub_log_sPTB<-sub_log_sPTB[order(sub_log_sPTB$AptName),]

df_heatmap<-data.frame(mPTB=sub_log_mPTB$logFC,
                       sPTB=sub_log_sPTB$logFC)

rownames(df_heatmap)=sub_log_mPTB$EntrezGeneID


heatmap_3<-heatmap_plot(df_heatmap)


# Save Plots----

## PLSDA Plots----
ggsave(sub_mPTB_sPTB_OE_neg_plot,
       path="../final_plots/figure3",
       width=7,height=8,units="in",
       filename="fig3C_Reactome_barplots.pdf")

## Heatmaps----
#Heatmap 1
pdf(file="../final_plots/figure3D_classI_pathway_heatmap.pdf",height=8,width=3)
draw(heatmap_1, heatmap_legend_side = "right", #annotation_legend_side = "bottom", 
     column_title="Proteins Associated with\nClass I MHC mediated antigen\nprocessing & presentation",
     column_title_gp=grid::gpar(fontface="bold"))
dev.off()

#Heatmap 2
pdf(file="../final_plots/figure3D_IL4andIL13_pathway_heatmap.pdf",height=8,width=3)
draw(heatmap_2, heatmap_legend_side = "right", #annotation_legend_side = "bottom", 
     column_title="Proteins Associated with\nInterleukin-4 and\nInterleukin-13 signaling",
     column_title_gp=grid::gpar(fontface="bold"))
dev.off()

#Heatmap 3
pdf(file="../final_plots/figure3D_IGF_pathway_heatmap.pdf",height=8,width=3)
draw(heatmap_3, heatmap_legend_side = "right", #annotation_legend_side = "bottom", 
     column_title="Proteins Associated with\nRegulation of IGF transport\nand uptake by IGFBPs",
     column_title_gp=grid::gpar(fontface="bold"))
dev.off()


# Save CSV/Excel Files ----


## All Pathway Results----
write_xlsx(sub_pathways_list,
           path="../../Supplemental_data_revisions/Sup_data1.xlsx")

## Main Figures----
save_csv_pathway = ALL_pathway_groups_merged[,c(2,9,10,11)]
rownames(save_csv_pathway) = NULL
write.csv(save_csv_pathway,"../../Supplemental_data_revisions/fig3b_pathway_barplot.csv")

heatmap_1_csv = heatmap_1@matrix[row_order(heatmap_1),]
write.csv(heatmap_1_csv,"../../Supplemental_data_revisions/fig3d_class_I_heatmap.csv")

heatmap_2_csv = heatmap_2@matrix[row_order(heatmap_2),]
write.csv(heatmap_2_csv,"../../Supplemental_data_revisions/fig3d_interleukin_heatmap.csv")

heatmap_3@csv = heatmap_3@matrix[row_order(heatmap_3),]
write.csv(heatmap_4_csv,"../../Supplemental_data_revisions/fig3d_igf_transport_heatmap.csv")


