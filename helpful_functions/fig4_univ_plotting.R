# this is going to be used both by Figure 4C and 4D Univ
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(rstatix)
  library(ggprism)
})

LRT_univ_fig3_stats<-function(apt_name, df_int){
  temp_feat = sym(apt_name)
  equ_comp = as.formula(paste0("`", apt_name,"`","~","Type"))
  stat_test <- df_int %>%
    group_by(bmi30) %>%
    wilcox_test(data=.,formula = equ_comp,paired=FALSE) %>%
    adjust_pvalue(method="bonferroni") %>%
    add_significance("p.adj")
  
  #removing 'ns' symbols
  stat_test <- remove_ns(stat_test) 
  
  #doing if statement since some features may not have any significant comparisons 
  if(nrow(stat_test)>0){
    stat_test <- stat_test %>%
      add_xy_position(x = "Type")
  }
  else{
    stat_test = NULL
  }
  
  return(stat_test)
}

LRT_univ_fig3_plot<-function(apt_name, 
                             df_int,
                             stat_test = NULL, #some features may not be included
                             apt_conv, #the t_test_v3 conv
                             color_groups) {
  
  temp_feat = sym(apt_name)
  
  #convert From Apt to Gene Symbol
  gene_symbol = subset(apt_conv, AptName %in% apt_name)$EntrezGeneSymbol
  
  temp_plot <- ggplot(df_int, aes(x=Type, y=!!temp_feat))+
    #plotting aesthetics
    geom_violin(draw_quantiles = c(0.25,0.5,0.75),aes(fill=Type),alpha = 0.5)+
    geom_point(aes(fill=Type),pch=21,size=1)+
    scale_fill_manual(values = color_groups)+
    facet_grid(~bmi30)+
    theme_classic()+
    #adding labels 
    ggtitle(gene_symbol,
            subtitle = "Stratified by BMI (kg/m2)")+
    ylab("log2(RFU)")+
    #fixing figure parameters
    theme(plot.title=element_text(hjust=0.5,face="bold",size=12),
          plot.subtitle = element_text(hjust=0.5,size=10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8),
          axis.text = element_text(size = 7),
          strip.text = element_text(face="bold",size = 9),
          #strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.3),
          legend.position = "none",
          axis.line = element_blank())
  
  
  #adding stats if applies
  if(!(is.null(stat_test))){
    temp_plot = temp_plot + add_pvalue(stat_test,
                                       label = "p.adj.signif",
                                       xmin = "xmin",
                                       bracket.size = 0.3,
                                       size = 1.25,
                                       #remove.bracket = TRUE,
                                       #size=1,
                                       xmax = "xmax",
                                       tip.length = 0)
  }
  
  return(temp_plot)
  
  
}