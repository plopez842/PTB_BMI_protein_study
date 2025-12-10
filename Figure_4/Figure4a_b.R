# Purpose of Script ---- 
# Figures for Figure 4A and B

# Preparing Packages and Functions----
set.seed(1010)
source("../helpful_functions/Systems_Serology_PAL.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(VennDiagram)
})

# Load Variables----
all_comp_limma<-readRDS(file="../output_RDS_objs/pairwise_comp_volc_plot_list_df.rds")

# Fig 3A - Volcano Plots----

comp_int <- c("g30Vsl30_ctrl",
              "g30Vsl30_mptb",
              "g30Vsl30_sptb")

new_title <- c( "Ctrl: >=30 vs <30 BMI",
                "mPTB: >=30 vs <30 BMI",
                "sPTB: >=30 vs <30 BMI")

sub_title<- c("Total DEP: 284",
              "Total DEP: 4",
              "Total DEP: 48")

volc_plot<- list()
myCol=c("#868686FF","#CD534CFF","#0073C2FF")

for (i in 1:length(comp_int)){
  title_name<-new_title[i]
  subtitle_name = sub_title[i]
  temp_comp_int <- comp_int[i]
  
  temp_col=myCol[i]
  
  temp_df<-all_comp_limma[[temp_comp_int]]
  
  temp_df<-temp_df[order(temp_df$phenotype,decreasing=TRUE),]
  
  sub_df<-subset(temp_df,phenotype %in% "yes")
  keep_labels<-sub_df[!(duplicated(sub_df$label)),]
  keep_labels<-keep_labels[order(abs(keep_labels$logFC),decreasing=TRUE),]
  keep_top_15<-rownames(keep_labels)[1:15]
  temp_df[!(rownames(temp_df) %in% keep_top_15),]$label=NA
  
  temp_plot<-ggplot(temp_df,aes(x=logFC,y=-log10(adj.P.Val),label=label))+
    geom_point(aes(color=phenotype,alpha=phenotype),size=0.75)+
    geom_text_repel(box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"),size=2.5,
                    max.overlaps = 20)+
    scale_color_manual(values=c("no"="gray","yes"=temp_col))+
    theme_classic()+
    
    #add in alpha
    scale_alpha_manual(values=c("no"=0.5,"yes"=1))+
    
    ylab("-log10(FDR)")+
    xlab("log2FC")+
    #
    xlim(-2,2)+
    ylim(0,2.5)+
    
    ggtitle(title_name,subtitle = subtitle_name)+
    geom_vline(xintercept=0.58,linetype="dashed",color="darkgray")+
    geom_vline(xintercept=-0.58,linetype="dashed",color="darkgray")+
    geom_hline(yintercept=0.602,linetype="dashed",color="darkgray")+
    theme(legend.position = "none",
          axis.title=element_text(face="bold",size=8),
          plot.title=element_text(hjust=0.5,face="bold",size=10),
          plot.subtitle = element_text(hjust=0.5,size=9.5),
          axis.text = element_text(size=6))
  volc_plot[[i]]<-temp_plot
}
names(volc_plot)<-comp_int


# Figure 3B - Venn Diagram----
sig_bmi_control<-subset(all_comp_limma$g30Vsl30_ctrl,phenotype %in% "yes")
sig_bmi_mptb<-subset(all_comp_limma$g30Vsl30_mptb,phenotype %in% "yes")
sig_bmi_sptb<-subset(all_comp_limma$g30Vsl30_sptb,phenotype %in% "yes")

## Plot and Save----
myCol=c("ctrl"="#868686FF","mptb"="#CD534CFF","sptb"="#0073C2FF")
venn.diagram(
  x=list(rownames(sig_bmi_control),rownames(sig_bmi_mptb),rownames(sig_bmi_sptb)),
  category.names = c("Ctrl","mPTB","sPTB"),
  filename="../final_plots/figure4/Fig4B_Venn_Diagram.png",
  output=TRUE,
  
  #main font
  main="Total Overall of # DEP in each\n>=30 vs <30 BMI Comparison",
  main.fontfamily = "sans",
  main.fontface = "bold",
  main.cex=0.6,
  main.pos = c(0.5, 0.9),
  imagetype="png" ,
  height = 480 , 
  width = 500 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  #fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.default.pos = "outer",
  cat.pos = c(-135, 135, 0),  # Adjust the category positions
  cat.dist = c(0.04, 0.035, 0.03),
  cat.fontfamily = "sans",
  rotation=180,
  scaled=FALSE,
  margin = 0.2
)


# Save Figures----
path_int<-"../final_plots/figure4"

lapply(names(volc_plot), function(x) ggsave(volc_plot[[x]],path=path_int,
                                            file=paste0("fig4A_",x,".pdf"),
                                            units="in",width=2.75,height=2.75))

# Save CSV files----
ctrl_only = all_comp_limma$g30Vsl30_ctrl

df_ctrl = data.frame(protein = ctrl_only$EntrezGeneID,
                     log2fc = ctrl_only$logFC,
                     FDR = ctrl_only$adj.P.Val)



mPTB_only = all_comp_limma$g30Vsl30_mptb
df_mptb = data.frame(protein = mPTB_only$EntrezGeneID,
                     log2fc = mPTB_only$logFC,
                     FDR = mPTB_only$adj.P.Val)

sPTB_only = all_comp_limma$g30Vsl30_sptb
df_sptb = data.frame(protein = sPTB_only$EntrezGeneID,
                     log2fc = sPTB_only$logFC,
                     FDR = sPTB_only$adj.P.Val)


write.csv(df_ctrl,"../../Supplemental_data_revisions/fig4a_ctrl_volc_plot.csv")
write.csv(df_mptb,"../../Supplemental_data_revisions/fig4a_mptb_volc_plot.csv")
write.csv(df_sptb,"../../Supplemental_data_revisions/fig4a_sptb_volc_plot.csv")

