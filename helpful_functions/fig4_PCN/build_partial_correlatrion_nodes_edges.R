

build_partial_correlation_network <- function(partial_corr_matrix,threshold=0) {
  # Convert partial correlation matrix to a standard matrix (if not already)
  partial_corr_matrix <- as.matrix(partial_corr_matrix)
  
  # Convert matrix to long format
  partial_corr_long <- melt(partial_corr_matrix)
  colnames(partial_corr_long) <- c("Node1", "Node2", "Correlation")
  
  # Remove diagonal elements and keep only one half of the matrix (for undirected graph)
  partial_corr_long <- partial_corr_long %>%
    filter(Node1 != Node2) %>%
    filter(as.numeric(Node1) > as.numeric(Node2))  # Keep only one half of the matrix
  
  # Set partial correlations less than 10^(-3) to zero
  partial_corr_long$Correlation[abs(partial_corr_long$Correlation) < 10^(-3)] <- 0
  
  partial_corr_long$Correlation[abs(partial_corr_long$Correlation) < threshold] <- 0
  
  # Remove edges where correlation is zero (no connection)
  partial_corr_long <- partial_corr_long %>%
    filter(Correlation != 0)
  
  
  
  # Create igraph object
  graph <- graph_from_data_frame(partial_corr_long, directed = FALSE)
  
  # Set edge weights based on the partial correlation values
  E(graph)$weight <- partial_corr_long$Correlation
  # Remove isolated nodes (nodes with no edges)
  graph <- delete.vertices(graph, V(graph)[degree(graph) == 0])
  return(graph)
}
