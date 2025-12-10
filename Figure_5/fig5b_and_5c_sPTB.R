# Purpsoe of Script----
## Prepare nodes and edges for cytoscape
## Most of the handling of aesthetics was done in Cytoscape, but need to properly setup nodes and edges files 

## Better to run this script in RStudio (and convert to RMarkdown) since it interacts with Cytoscape
# Load Libraries----

# Load Functions and Libraries----

functions_path_dir="../helpful_functions/fig4_PCN"
r_files <- list.files(path = functions_path_dir, pattern = "\\.R$", full.names = TRUE)
sapply(r_files,source)

suppressPackageStartupMessages({
  library(visNetwork)
  library(reshape2)
  library(RCy3)
})

# Load Variables----

full_corr = readRDS(file = "../output_RDS_objs/Fig5_PCN_int_objects/Fig5_sPTB_sampleCOV.RDS")
partial_corr = readRDS(file = "../output_RDS_objs/Fig5_PCN_int_objects/Fig5_sPTB_PCN.RDS")
protein_test_info = readRDS(file =  "../output_RDS_objs/fig5_protein_info.rds")

# Setup Nodes/Edges DF for Full_corr----
threshold <- 0.15

network_data_FULL <- build_partial_correlation_network(full_corr, threshold)
network_data_FULL<- toVisNetworkData(network_data_FULL)

# Extract nodes and edges
vis_DF_nodes_FULL <-network_data_FULL$nodes
vis_DF_edges_FULL <- network_data_FULL$edges

df_merge_vis_nodes_FULL = vis_DF_nodes_FULL %>%
  left_join(protein_test_info,by=c("label"="node_name"))

## Adjust Some Aesthetics
df_merge_vis_nodes_FULL$node_transparency<-1
df_merge_vis_nodes_FULL[df_merge_vis_nodes_FULL$protein_or_covariate=="covariate",]$node_transparency = 255

df_merge_vis_nodes_FULL$node_height<-35
df_merge_vis_nodes_FULL[df_merge_vis_nodes_FULL$protein_or_covariate=="covariate",]$node_height = 200


# Run Cytoscape Directly from the code below----
nodes_edges_FULL<-data.frame(source = vis_DF_edges_FULL$from,
                             target = vis_DF_edges_FULL$to,
                             weight = vis_DF_edges_FULL$weight)

createNetworkFromDataFrames(nodes=df_final_FULL ,
                            edges = nodes_edges_FULL)

# Setup Nodes/Edges DF for Partial Corr----

network_data_partial <- build_partial_correlation_network(partial_corr, threshold)
network_data_partial<- toVisNetworkData(network_data_partial)

# Extract nodes and edges
vis_DF_nodes_PARTIAL <-network_data_partial$nodes
vis_DF_edges_PARTIAL <- network_data_partial$edges

df_merge_vis_nodes_PARTIAL = vis_DF_nodes_PARTIAL %>%
  left_join(protein_test_info,by=c("label"="node_name"))

## Adjust Some Aesthetics
df_merge_vis_nodes_PARTIAL$node_transparency<-1
df_merge_vis_nodes_PARTIAL[df_merge_vis_nodes_PARTIAL$protein_or_covariate=="covariate",]$node_transparency = 255

df_merge_vis_nodes_PARTIAL$node_height<-35
df_merge_vis_nodes_PARTIAL[df_merge_vis_nodes_PARTIAL$protein_or_covariate=="covariate",]$node_height = 200


# Run Cytoscape Directly from the code below----
nodes_edges_PARTIAL<-data.frame(source = vis_DF_edges_PARTIAL$from,
                             target = vis_DF_edges_PARTIAL$to,
                             weight = vis_DF_edges_PARTIAL$weight)

createNetworkFromDataFrames(nodes=df_final_PARTIAL ,
                            edges = nodes_edges_PARTIAL)