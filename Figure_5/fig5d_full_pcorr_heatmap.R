# Purpose of Script----
# Construct Figure 5D Heatmap

# Load Libraries----
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(colorRamp2)
})

# Load Variables----

saveRDS(full_corr_GENES,"../output_RDS_objs/Fig5_PCN_int_objects/Fig5_mPTB_sampleCOV.RDS")
saveRDS(partial_corr_GENES,"../output_RDS_objs/Fig5_PCN_int_objects/Fig5_mPTB_PCN.RDS")
mPTB_full_corr<-readRDS(file="../output_RDS_objs/Fig5_PCN_int_objects/Fig5_mPTB_sampleCOV.RDS")
sPTB_full_corr<-readRDS(file="../output_RDS_objs/Fig5_PCN_int_objects/Fig5_sPTB_sampleCOV.RDS")

mPTB_PCN<-readRDS(file="../output_RDS_objs/Fig5_PCN_int_objects/Fig5_mPTB_PCN.RDS")
sPTB_PCN<-readRDS(file="../output_RDS_objs/Fig5_PCN_int_objects/Fig5_sPTB_PCN.RDS")

# Work with Partial Corr First


PCN_corr<-matrix(nrow=149,ncol=4)
rownames(PCN_corr) = rownames(mPTB_full_corr)[3:151]
colnames(PCN_corr) = c("dirr_with_BMI_full_mPTB",
                       "dirr_with_BMI_full_sPTB",
                       "dirr_with_PTB_full_mPTB",
                       "dirr_with_PTB_full_sPTB")

PCN_corr[,1] = mPTB_PCN[3:151,1]
PCN_corr[,2] = sPTB_PCN[3:151,1]
PCN_corr[,3] = mPTB_PCN[3:151,2]
PCN_corr[,4] = sPTB_PCN[3:151,2]

## Subset to Proteins that have at least one nonzro Partial Corr to either PTB or BMi
which_proteins_to_keep = rownames(PCN_corr[!apply(PCN_corr,1,function(x) all(is.na(x))),])
## Prepare Annotations
group_annotation_Partial = HeatmapAnnotation("Network" = c("mPTB vs BMI","sPTB vs BMI","mPTB vs BMI","sPTB vs BMI"),
                                             col = list(Network = c("mPTB vs BMI" ="#CD534CFF" ,"sPTB vs BMI" = "#0073C2FF")),
                                             show_annotation_name = FALSE,
                                             show_legend = TRUE,
                                             border=TRUE)

corr_col = colorRamp2(c(-1,0,1),c("#1207A3","#FFFEFE","#BB0103"))

## Plot
sub_heatmap_PCN = PCN_corr[which_proteins_to_keep,]
sub_PCN<-Heatmap(sub_heatmap_PCN ,column_split=factor(rep(c("Partial Corr\nwith BMI","Partial Corr\nwith PTB"),each=2)),
                 column_title_gp = gpar(fontsize = 8),
                 show_row_names = FALSE,
                 top_annotation = group_annotation_Partial,
                 name = "Partial Corr",
                 col=corr_col,border_gp = gpar(col = "black"),
                 cluster_rows = FALSE,cluster_columns = FALSE,
                 heatmap_legend_param = list(direction = "horizontal",legend_side="bottom",title_position = "topcenter",legend_width = unit(20,"mm"),legend_height = unit(1,"mm"),at = c(-1,0,1)),
                 row_names_side = "left",
                 width = ncol(full_corr)*unit(7 ,"mm"),
                 height = nrow(full_corr)*unit(2.5,"mm"),
                 show_column_names = FALSE,
                 row_names_gp = gpar(fontsize = 6))


# Work with full Corr ----
full_corr<-matrix(nrow=149,ncol=4)
rownames(full_corr) = rownames(mPTB_full_corr)[3:151]
colnames(full_corr) = c("dirr_with_BMI_full_mPTB",
                        "dirr_with_BMI_full_sPTB",
                        "dirr_with_PTB_full_mPTB",
                        "dirr_with_PTB_full_sPTB")

full_corr[,1] = mPTB_full_corr[3:151,1]
full_corr[,2] = sPTB_full_corr[3:151,1]
full_corr[,3] = mPTB_full_corr[3:151,2]
full_corr[,4] = sPTB_full_corr[3:151,2]

## Prepare Annotations----
group_annotation_FULL = HeatmapAnnotation("Network" = c("mPTB vs BMI","sPTB vs BMI","mPTB vs BMI","sPTB vs BMI"),
                                          col = list(Network = c("mPTB vs BMI" ="#CD534CFF" ,"sPTB vs BMI" = "#0073C2FF")),
                                          show_annotation_name = FALSE,
                                          show_legend = FALSE,
                                          border=TRUE)

corr_col = colorRamp2(c(-1,0,1),c("#1207A3","#FFFEFE","#BB0103"))

## Plot only Full Heatmap----
sub_heatmap_full = full_corr[which_proteins_to_keep,]
sub_full<-Heatmap(sub_heatmap_full,column_split=factor(rep(c("Corr with\nBMI","Corr with\nwith PTB"),each=2)),
                  column_title_gp = gpar(fontsize = 8),
                  top_annotation = group_annotation_FULL,
                  name = "Full Corr",
                  col=corr_col,border_gp = gpar(col = "black"),
                  cluster_rows = FALSE,cluster_columns = FALSE,
                  heatmap_legend_param = list(direction = "horizontal",legend_side="bottom",title_position = "topcenter",legend_width = unit(20,"mm"),legend_height = unit(1,"mm"),at = c(-1,0,1)),
                  row_names_side = "left",
                  width = ncol(sub_heatmap_full)*unit(7 ,"mm"),
                  height = nrow(sub_heatmap_full)*unit(2.5,"mm"),
                  #legend_width = unit(1,"cm"),
                  show_column_names = FALSE,
                  #adjust the param
                  row_names_gp = gpar(fontsize = 6))

# Save Plots----
sub_merged<- sub_full+sub_PCN
pdf(file="../output_plots/figure5/Fig5D_Full_Partial_HEATMAP_merged.pdf",height=9,width=6)
draw(sub_merged,heatmap_legend_side="bottom",column_title = "Correlation vs Partial correlation of Proteins\nto BMI and PTB status")
dev.off()

# Save CSV----
write.csv(sub_heatmap_full,"../../Supplemental_data_revisions/fig5b_full_heatmap.csv")
write.csv(sub_heatmap_PCN,"../../Supplemental_data_revisions/fig5b_partial_heatmap.csv")


