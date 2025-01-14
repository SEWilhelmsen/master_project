# Run create network plot 
# Silje Wilhelmsen


# Source necessary
source("C:/Users/siljeew/snRNAseq/Visualization/load_libraries_run_create_network_plot.R")
source("C:/Users/siljeew/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/snRNAseq/Visualization/create_network_plot_functions.R")

# Set parameters
#############################################################################
output_dir_network <- "C:/Users/siljeew/snRNAseq/Plots/Network"
# Ensure the directory exists
if (!dir.exists(output_dir_network)) {
  dir.create(output_dir_network, recursive = TRUE)
}

# Define specific processes for your data set
go_processes_of_interest <- c('glycolysis')

# Define the specific processes you want to include as nodes
specific_process_names <- c("glycolysis", "pdh", "etc", "fao", "tca", "ketone_catabolism", "protein_transmembrane_transport", "chylomicron_clearance", "atp_synthesis", "uncoupling_proteins", "ketone_metabolic_process")

file_path <- "C:/Users/siljeew/snRNAseq/tmp/mouse_6h_vcm.Rds"

threshold <- 0.7




# Load data and subset based on threshold
##############################################################################
# Load necessary data and prepare the subset
subset_expression_matrix <- load_and_preprocess_data(file_path, go_processes_of_interest, define_go_process, get_genes_for_go_process)

# Compute correlation matrix
correlation_matrix <- compute_correlation_matrix(subset_expression_matrix)
genes <- rownames(correlation_matrix)

# Create adjacency matrix
adjacency_matrix <- create_adjacency_matrix(correlation_matrix, genes, threshold)



# Network construction
################################################################################
# Retrieve GO terms for genes in your network
gene_symbols <- rownames(subset_expression_matrix)  # Use rownames of filtered expression data
gene_go_pairs <- get_go_terms_for_genes(gene_symbols) 
process_nodes <- names(define_go_process)
gene_process_edges <- data.frame(from = gene_go_pairs$GOALL, to = gene_go_pairs$SYMBOL)



# Filter and subset nodes
#############################################################################
process_nodes <- names(define_go_process)
gene_process_edges <- data.frame(from = gene_go_pairs$GOALL, to = gene_go_pairs$SYMBOL)

# Update edge types
node_labels <- setNames(c(process_nodes), define_go_process)

# Confirm all required vertices are in your edge frame
all_vertices <- unique(c(gene_process_edges$from, gene_process_edges$to))
vertices_needed <- unique(c(process_nodes, genes))

# Identify any missing vertex labels
missing_vertices <- setdiff(all_vertices, vertices_needed)
print(missing_vertices)  # This should show any vertices missing

# Verify process nodes and genes consistency
vertices_corrected <- unique(c(process_nodes, genes, missing_vertices))


# Inspect number of unique processes
unique_processes <- unique(gene_process_edges$from) 
number_of_unique_processes <- length(unique_processes)  
cat("Number of unique processes:", number_of_unique_processes, "\n")

# Extract the corresponding GO IDs for these processes
included_go_process_ids <- unlist(define_go_process[specific_process_names])
print(included_go_process_ids)

# Filter the gene_process_edges data frame
filtered_gene_process_edges <- gene_process_edges[gene_process_edges$from %in% included_go_process_ids, ]

# Print number of unique filtered processes
unique_filtered_processes <- unique(filtered_gene_process_edges$from)
number_of_filtered_processes <- length(unique_filtered_processes)
cat("Number of unique filtered processes:", number_of_filtered_processes, "\n")

# Inspect filtered data frame
print(head(filtered_gene_process_edges))


# Create network graph
##############################################################################

# Create the undirected graph
graph <- graph_from_data_frame(filtered_gene_process_edges, directed = FALSE, 
                               vertices = unique(c(filtered_gene_process_edges$from, filtered_gene_process_edges$to)))

# Set layout of edges/lines
E(graph)$color <- "grey"
E(graph)$arrow.size <- 0.5


# Set process names 
process_names <- setNames(names(define_go_process), unlist(define_go_process))
print(process_names)

# Set vertex labels for clarity
V(graph)$label <- ifelse(V(graph)$name %in% filtered_gene_process_edges$to,
                         V(graph)$name,  # Gene names as labels
                         process_names[V(graph)$name])  # Process names as labels

print(head(V(graph)$label))

# Determine vertex type for styling
V(graph)$type <- ifelse(V(graph)$name %in% filtered_gene_process_edges$to, "gene", "process")

# Style vertices: Set gene nodes to be transparent and size zero
V(graph)$color <- ifelse(V(graph)$type == "process", "bisque1", "white")  
V(graph)$size <- ifelse(V(graph)$type == "process", 10, 7) 

# Adjust label size for better visibility if needed
V(graph)$label.cex <- ifelse(V(graph)$type == "process", 0.8, 0.7)

# Select an appropriate layout
layout <- layout_with_fr(graph)

# Plot the graph
plot(graph, 
     layout = layout, 
     vertex.label.cex = .6, 
     vertex.frame.color = NA, 
     main = "Gene-Process Network", 
     vertex.label.color = "black", 
     edge.curved = .1, 
     edge.arrow.size = .3, 
     edge.width = .7)

