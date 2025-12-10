
suppressPackageStartupMessages({
  library(DOSE)
  
})
extract_geneSets <- function(x, n) {
  n <- update_n(x, n)
  geneSets <- geneInCategory(x) ## use core gene for gsea result
  y <- as.data.frame(x)
  geneSets <- geneSets[y$ID]
  names(geneSets) <- y$Description
  if (is.numeric(n)) {
    return(geneSets[1:n])
  }
  return(geneSets[n]) ## if n is a vector of Description
  
  
}


update_n <- function(x, showCategory) {
  if (!is.numeric(showCategory)) {
    return(showCategory)
  }
  
  ## geneSets <- geneInCategory(x) ## use core gene for gsea result
  n <- showCategory
  if (nrow(x) < n) {
    n <- nrow(x)
  }
  
  return(n)
}


apt_from_pathways<-function(pathway_int,
                            merged_pathway_groups,
                            corr_vals_list,
                            which_group){    
  
  acceptable_groups<-c("mPTB_pos","mPTB_neg","sPTB_neg")
  
  if(which_group %in% acceptable_groups){
    protein_heatmap_start1<-setReadable(merged_pathway_groups[[which_group]], 'org.Hs.eg.db', 'UNIPROT')
    protein_heatmap_start1@pvalueCutoff=1
    # protein_heatmap_start1@qvalueCutoff=1
    geneSets1<-extract_geneSets(protein_heatmap_start1,nrow(protein_heatmap_start1@result))
    gene_set_vect1<-geneSets1[[pathway_int]]
    if(is.null(gene_set_vect1)){
      print("pathway name bad")
    }
    else{
      df_apt_int<-subset(corr_vals_list[[which_group]],Single_gene %in% gene_set_vect1)
      df_apt_int_unique <- df_apt_int %>%
        group_by(genes) %>%
        slice_max(order_by = fscore, n = 1) %>%
        ungroup()
      
      
      apt_vect<-df_apt_int_unique$protein
      
      #write a case for duplication
      
      
      return(apt_vect)
    }
    
  }
  else{
    print("Comparison group not in list")
  }
  
}

heatmap_plot<-function(df_heatmap){
  df_heatmap<-as.matrix(df_heatmap)
  row_ha = HeatmapAnnotation(#BMI = sub_merged_data$bmi30,
    Status=c("mPTB","sPTB"),
    col = list(Status=c("mPTB"="#CD534CFF","sPTB"="#0073C2FF","ctrl"="#868686FF")),
    border=TRUE,show_legend=FALSE,
    show_annotation_name = FALSE)
  #log2_col= colorRamp2(c(min(df_heatmap), 0, max(df_heatmap)), c("#1207A3", "#FFFEFE", "#BB0103"))
  log2_col= colorRamp2(c(-0.5, 0, 0.5), c("#1207A3", "#FFFEFE", "#BB0103"))
  log2_class1<-Heatmap(df_heatmap,
                       show_heatmap_legend = FALSE,
                       # column_title= title,
                       #column_title_gp = gpar(fontface = "bold"),
                       name="log2FC\nvs Ctrl",
                       border_gp=gpar(col="black",lwd=1),
                       heatmap_legend_param = list(legend_height = unit(6, "cm"),border = "black"),
                       top_annotation = row_ha,
                       col=log2_col,
                       show_row_dend = FALSE,
                       cluster_columns = FALSE,
                       show_column_dend = FALSE,
                       width = ncol(df_heatmap)*unit(8 ,"mm"),
                       height = nrow(df_heatmap)*unit(6,"mm"),
                       column_names_side="top",
                       row_names_side="left"
  )
  
  return(log2_class1)
  
}