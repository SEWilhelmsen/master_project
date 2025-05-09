# Create network plot for genes and process
# Silje Wilhelmsen

library(igraph)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(GO.db)
library(dplyr)


# Source functions
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")

# Set parameters
##############################################################################
output_dir_network <- "C:/Users/siljeew/Master_project/snRNAseq/Plots/Network"

# Ensure the directory exists
if (!dir.exists(output_dir_network)) {
  dir.create(output_dir_network, recursive = TRUE)
}

# Define main process of interest and additional processes and the time point of interest
go_processes_of_interest <- c('glycolysis')
specific_process_names <- c("glycolysis", "pdh", "b_oxidation", "tca", "ketone_catabolism")
time_point_of_interest <- "6 Hours"



# Load and prepare data 
###################################################
# Data frame containing gene expression 
mouse_vcm_all_genes_avg_logfc <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/Data/mouse_vcm_all_genes_avg_logFC.rds")

# Change column names to better readability
mouse_vcm_all_genes_avg_logfc <- mouse_vcm_all_genes_avg_logfc %>%
  rename("6 Hours" = logFC6_hours, 
         "12 Hours" = logFC12_hours, 
         "1 Day" = logFC1_day, 
         "3 Days" = logFC3_days, 
         "1 Week" = logFC1_week, 
         "3 Weeks" = logFC3_weeks)

# Subset for time point 
data_filtered <- mouse_vcm_all_genes_avg_logfc[, time_point_of_interest, drop = FALSE]

# Set row names (gene names) if not automatically retained
rownames(data_filtered) <- rownames(mouse_vcm_all_genes_avg_logfc)

# Remove rows with NA values
data_filtered <- data_filtered[!is.na(data_filtered[, time_point_of_interest]), , drop = FALSE]



# Prepare filtering 
####################################################
# Define GO processes and filter genes
genes_of_interest_vector <- get_genes_of_interest(go_processes_of_interest, define_go_process, get_genes_for_go_process)
genes_of_interest_unlisted <- unlist(genes_of_interest_vector, use.names = TRUE)

# Keep the rownames that exist in both data_filtered and genes_of_interest
filtered_rows <- rownames(data_filtered) %in% genes_of_interest_unlisted

# Keep a data frame with only common rows from filtered_rows
data_filtered <- data_filtered[filtered_rows, , drop = FALSE]
head(data_filtered)
nrow(data_filtered)


# Prepare association data frame
##########################################################################
# Get unique list of genes across all processes
genes <- unique(unlist(genes_of_interest_vector))

# Function to get all GO terms for a list of genes
get_go_terms_for_genes <- function(genes) {
  go_terms <- AnnotationDbi::select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", columns = c("GOALL", "SYMBOL"))
  go_terms <- go_terms[!is.na(go_terms$GOALL), ]  # Drop entries without GO term
  unique(go_terms)
}

# Retrieve GO terms for genes of interest
gene_go_pairs <- get_go_terms_for_genes(genes)
nrow(gene_go_pairs) # Number of rows: 40798

# Create edges with genes to GO terms
edges <- data.frame(from = gene_go_pairs$SYMBOL, to = gene_go_pairs$GOALL)
nrow(edges) # Prints 40798

# Collect the specific_process_names from the define_go_process list
go_process_for_nodes <- unlist(define_go_process[specific_process_names])

# Filter the edges to only include connections to the defined processes
filtered_edges <- edges[edges$to %in% go_process_for_nodes, ]
print(nrow(filtered_edges))

# Get all unique vertices needed for this subset (genes plus chosen processes)
filtered_vertex_names <- unique(c(filtered_edges$from, filtered_edges$to))
vertices <- data.frame(name = filtered_vertex_names)
head(vertices)

# Make sure only the genes from data_filtered are present in vertices
genes_with_expression <- rownames(data_filtered)
genes_in_vertices <- vertices$name[!(grepl("^GO:", vertices$name))]
go_terms_in_vertices <- vertices$name[grepl("^GO:", vertices$name)]

# Filter vertices to only include the genes and GO terms
filtered_vertex_names <- unique(c(genes_with_expression, go_terms_in_vertices))

# Create the filtered vertices data frame
vertices <- data.frame(name = filtered_vertex_names)
head(vertices)

# Inspect results
print(nrow(vertices))
print(head(vertices))



# Create nodes of genes and processes
############################################################################
# Define GO IDs and retrieve their names from GO.db
go_ids <- unlist(define_go_process)
go_terms <- AnnotationDbi::select(GO.db,
                                  keys = go_ids,
                                  keytype = "GOID",
                                  columns = c("GOID", "TERM"))
go_names <- setNames(go_terms$TERM, go_terms$GOID)
process_nodes <- names(define_go_process)
gene_process_edges <- data.frame(from = gene_go_pairs$GOALL, to = gene_go_pairs$SYMBOL)


# Filter processes
included_go_process_ids <- unlist(define_go_process[specific_process_names])
filtered_gene_process_edges <- gene_process_edges[gene_process_edges$from %in% included_go_process_ids, ]

# Output a summary of filtered edges
unique_filtered_processes <- unique(filtered_gene_process_edges$from)
number_of_filtered_processes <- length(unique_filtered_processes)
cat("Number of unique filtered processes:", number_of_filtered_processes, "\n") # Prints 6
print(head(filtered_gene_process_edges))

# Extract the 'name' column from 'vertices', which contains all valid vertex names
valid_gene_vertices <- vertices$name

# Filter 'filtered_gene_process_edges' to keep only those rows with 'to' values in 'valid_gene_vertices'
filtered_gene_process_edges <- filtered_gene_process_edges[
  filtered_gene_process_edges$to %in% valid_gene_vertices, 
]

# Verify the updated edges
print(head(filtered_gene_process_edges))
nrow(filtered_gene_process_edges)  # Confirm the number of edges after filtering



# Create and save plot
############################################################################
# Create graph using the node data
graph <- graph_from_data_frame(filtered_gene_process_edges, directed = FALSE, vertices = vertices)
# Identify type: process or gene
V(graph)$type <- ifelse(V(graph)$name %in% names(go_names), "process", "gene")


# Label vertices
V(graph)$label <- sapply(seq_along(V(graph)), function(i) {
  vertex_name <- V(graph)$name[i]
  if (V(graph)$type[i] == "process") {
    process_label <- go_names[vertex_name]
    return(ifelse(!is.na(process_label) && process_label != "", process_label, vertex_name))
  } else {
    return(vertex_name)  # For genes, show gene name directly
  }
})

# Debugging outputs
print(table(V(graph)$type))
print(head(V(graph)$label))

# Visual styling for nodes
V(graph)$color <- ifelse(V(graph)$type == "process", "lightblue", "lightgreen")  #Colors of process and gene
V(graph)$size <- ifelse(V(graph)$type == "process", 8, 5)

# Use layout for graph
layout <- layout_with_fr(graph)


# Plot (basic igraph plot)
plot(graph, 
     layout = layout, 
     vertex.label.cex = 0.8, 
     vertex.frame.color = NA, 
     main = "Gene-Process Network Visualization", 
     vertex.label.color = "black", 
     edge.color = "gray",
     edge.curved = 0.1, 
     edge.arrow.size = 0.5, 
     edge.width = 0.7)


# Enhanced Plot with ggraph
create_network <- function(data_filtered, graph) {
  vertex_names <- V(graph)$name
  
  # Coloring vertices based on expression
  expression_colors <- sapply(V(graph)$name, function(name) {
    if (name %in% rownames(data_filtered)) {
      fold_change_value <- data_filtered[name, , drop = TRUE]
      if (!is.na(fold_change_value) && fold_change_value > 0) {
        return("rosybrown1")
      } else {
        return("lightblue1")
      }
    } else {
      return("grey")
    }
  })
  
  # Assign these colors to your vertices
  V(graph)$color <- expression_colors
  
  
  # Create advanced network plot using ggraph
  network_plot <- ggraph(graph, layout = 'fr') +
    geom_edge_link(aes(edge_alpha = 0.5), color = "grey", show.legend = FALSE) +
    geom_node_label(aes(label = label, fill = I(color)), fontface = "plain", color = "black", size = 6, label.size = 0, repel = TRUE) +
    scale_fill_identity() +
    theme_void() +
    labs(title = "Gene-Process Network of Differentially Expressed Genes")
  
  return(network_plot)
}


# Create plot with coloring based on expression data in the advanced layout
networkplot <- create_network(data_filtered, graph)
print(networkplot)



# Save plot 
################################################################################
file_name <- paste0(time_point_of_interest, "_network.png")
file_path <- file.path(output_dir_network, file_name)

# Save the plot to a file
png(file_path, width = 1000, height = 1000)
print(networkplot)
dev.off()

# Confirm the file path if needed
print(paste("Plot saved at:", file_path))
