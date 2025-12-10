
suppressPackageStartupMessages({
  library(dplyr) # for %>%
  library(rstatix) # for dunn_test 
  library(ggprism) # for adding pvalues to ggplot (add_pvalue)
})

# Fig 2d plots 
plotting_3way_PLSDA_univariates<-function(df, #unscaled dataframe
                                          feat, # protein of interest
                                          stat="Y", #whether or not to perform stat test
                                          ref_protein_Apt_mat, #protein Apt Info
                                          group_colors){ # colors of model 
  
  #initialize variables and plot
  protein_int=sym(feat)
  temp_plot<-ggplot(df,aes(x=Type,y=!!protein_int))+
    geom_violin(aes(fill=Type),draw_quantiles=c(0.25,0.5,0.75),alpha=0.3)+
    geom_point(aes(fill=Type),pch=21)+
    scale_fill_manual(guide="none",values=group_colors)+
    ggtitle(subset(ref_protein_Apt_mat,AptName %in% feat)$EntrezGeneSymbol)+
    theme_classic()+
    ylab("log2(RFU)")+
    theme(plot.title=element_text(hjust=0.5,face="bold"),
          axis.title.x=element_blank())
  
  #print(temp_plot)
  if(stat=="Y"){
    sub_eq<-as.formula(paste0("`", feat,"`","~","Type"))
    stat_test<-df %>%
      dunn_test(data=.,sub_eq,p.adjust.method="BH") 
    
    stat_test<-remove_ns(stat_test)
    #print(stat_test)
    if (nrow(stat_test)>0){
      print(stat_test)
      stat_test<-stat_test %>% add_xy_position(x="Type")
      

      temp_plot<-temp_plot+add_pvalue(stat_test,label="p.adj.signif", xmin="xmin", xmax="xmax",tip.length=0,
                                      label.size=4,step.increase=0.001)
    }
    else{
      temp_plot<-temp_plot
    }
    
  }
  return(temp_plot)
}