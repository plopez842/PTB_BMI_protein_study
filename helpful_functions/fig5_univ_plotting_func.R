
suppressPackageStartupMessages({
  library(ggplot2)
})

plotting_univariates<-function(df_of_data,
                               PCN_data,
                               seq_to_gene_df,
                               gene_int,
                               PTB_group = "mPTB"){
  # gene_int <- "UBE2G2"
  df_int = df_of_data
  pcn_itn = PCN_data
  
  
  #filtered_sub_t_tests
  color_and_fill_vals <- c("mPTB"="#CD534CFF","Ctrl"="#868686FF","sPTB"="#0073C2FF")
  seq_of_int<-subset(seq_to_gene_df, Single_gene %in% gene_int )$AptName
  df_int[df_int$Type2==1,]$Type2 =PTB_group
  df_int[df_int$Type2==0,]$Type2 = "Ctrl"
  
  print(gene_int)
  #prep text for plot
  pcorr_bmi = pcn_itn['BMI',gene_int]
  pcorr_bmi = format(pcorr_bmi , digits = 2)
  
  pcorr_ptb = pcn_itn['PTB_status',gene_int]
  pcorr_ptb = format(pcorr_ptb , digits = 2)
  
  pcorr_merged = paste0("p corr w/ BMI = ",pcorr_bmi,"\np.corr w/ PTB = ",pcorr_ptb)
  
  
  symb_of_int = sym(seq_of_int)
  
  # prepare the annotation
  # 1. Fit the linear model
  form <- as.formula(paste(seq_of_int, "~ bmi1"))
  model <- lm(form,data = df_int)
  
  # 2. Extract the slope
  slope <- coef(model)[["bmi1"]]
  
  # 3. Set x and hjust based on slope direction
  x_pos <- ifelse(slope > 0, 19, 51)  # left or right
  hjust_val <- ifelse(slope > 0, 0, 1)  # left-align or right-align
  
  # 4. Set y at the top of the y-axis range
  #y_max <- max(df_int[[as_string(symb_of_int)]], na.rm = TRUE)
  
  plot<- ggplot(df_int,aes(x=bmi1,y=!!symb_of_int,color=Type2))+
    geom_point(aes(fill=Type2),pch=21,alpha=0.2,size=0.75)+
    scale_fill_manual(values=color_and_fill_vals)+
    geom_smooth(aes(fill=Type2),method='lm', formula= y~x,alpha=0.3,level=0.7)+
    scale_color_manual(values=color_and_fill_vals)+
    xlim(18,52)+
    xlab("BMI (kg/m3)")+
    ylab(paste0("log2(RFU)"))+
    theme_classic()+
    theme(axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title = element_text(face="bold",size=10),
          axis.text = element_text(size=6),
          legend.position = "none")+
    annotate("text", x = x_pos, y = max(df_int[seq_of_int], na.rm = TRUE) - 0.1,
             label = pcorr_merged , hjust = hjust_val, size = 2)
  
  return(plot)
}