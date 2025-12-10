# Purpose of Script ---- 
# Figures for Figure 1A

# Loading Packages ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
})

# Laoding Objects ----
# these objects were prepared externally with clinical collaborators 
patient_info_levels = readRDS("../output_RDS_objs/merged_metadata_and_protein_V1.rds")

# Formating Obj for Diagram ----
diagram_obj<-patient_info_levels %>% 
  pivot_longer(
    cols= c("ga_exact","ga_blooddraw"),
    names_to=c("GA_type"),
    values_to=c("val")
    
  )
diagram_obj<-diagram_obj %>%
  group_by(Type) %>%
  arrange(val)

diagram_obj$studyid<-factor(diagram_obj$studyid,levels=unique(diagram_obj$studyid))
diagram_obj$GA_type<-factor(diagram_obj$GA_type,levels=c("ga_blooddraw","ga_exact"))
diagram_obj$GA_type<-factor(diagram_obj$GA_type,labels=c("At Sampling","At Delivery"))

diagram_obj$Type<-factor(diagram_obj$Type,levels=c("ctrl","meta_sptb","mtpb"))
diagram_obj$Type<-factor(diagram_obj$Type,labels=c("Control\nn=40","sPTB\nn=30","mPTB\nn=30"))

# Plotting obj ----
group_cols_diag<-c("Control\nn=40"="#868686FF","sPTB\nn=30"="#0073C2FF","mPTB\nn=30"="#CD534CFF")#0073C2FF


fig1A_plot<-ggplot(diagram_obj,aes(x=val,y=(studyid),group=studyid))+
  geom_point(aes(shape=GA_type,fill=Type))+
  scale_shape_manual(values = c(21,25)) +
  scale_fill_manual(guide="none",values=group_cols_diag)+
  geom_line()+
  theme_classic()+
  facet_wrap(vars(Type),
             ncol=1,
             strip.position="left",
             scales="free_y")+
  xlab("Gestation Age (weeks)")+
  geom_vline(xintercept = 28,linetype="dashed",colour="darkgrey")+
  theme(axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        legend.box.background = element_rect(colour = "black",size=0.5),
        legend.position = c(0.875, 0.1),
        strip.background = element_blank(),
        legend.title=element_blank(),
        strip.text.y=element_text(angle=0,face="bold",size=10),
        legend.text = element_text(size=8),
        axis.title.x=element_text(size=12),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

# Save Plot ----
ggsave(filename = "fig1a_diagram_of_Ns.pdf",temp_diag,units="in",
       path="../final_plots/figure1",
       width=5.5,
       height=6.5)