## Create heatmap from combined_avg_log2fc
# Silje Wilhelmsen

# source("C:/Users/siljeew/snRNAseq/Visualization/Uferdig/load_libraries_run_create_heatmap_combined_avg_log2fc.R")
# source("C:/Users/siljeew/snRNAseq/Analysis/define_go_terms_and_find_genes.R")


# # define_go_terms_and_find_genes.R
# ##############################################################################
# # Define GO processes
# define_go_process <- list(
#   glycolysis = "GO:0006096",
#   pdh = "GO:0045254",
#   fao = "GO:0019395",
#   tca = "GO:0006099",
#   etc = "GO:0022900",
#   chylomicron_clearance = "GO:0034382",
#   atp_synthesis = "GO:0042775",
#   etc_succinate_to_ubiquinone = "GO:0006121",
#   uncoupling_proteins = "GO:0017077",
#   fatty_acid_transport = "GO:0015908",
#   b_oxidation = "GO:0006635",
#   monosaccharide_transmembrane_transport = "GO:1905950",
#   protein_transmembrane_transport = "GO:0071806",
#   dna_binding_transcription_factor_activity = "GO:0003700",
#   translation_regulatory_activity = "GO:0045182",
#   ketone_metabolic_process = "GO:0042180",
#   hydroxybutyrate_dehydrogenase_activity = "GO:0003858",
#   regulation_of_ketone_metabolic_process = "GO:0010565",
#   ketone_catabolism = "GO:0046952"
# )
# 
# # Function to retrieve genes for a GO term
# get_genes_for_go_process <- function(go_term) {
#   entrez_ids <- AnnotationDbi::select(org.Hs.eg.db,
#                                       keys = go_term,
#                                       keytype = "GOALL",
#                                       columns = "ENTREZID")[, "ENTREZID"]
# 
#   gene_symbols <- AnnotationDbi::select(org.Hs.eg.db,
#                                         keys = entrez_ids,
#                                         keytype = "ENTREZID",
#                                         columns = "SYMBOL")[, "SYMBOL"]
# 
#   # Normalize gene symbols
#   gene_symbols_upper <- toupper(gene_symbols)
# 
#   # Return unique gene symbols related to the GO term
#   unique(gene_symbols_upper)
# }
# 
# # Function to retrieve genes organized by process
# get_genes_of_interest <- function(go_terms, define_go_process, get_genes_for_go_process) {
#   # Initialize a named list to store genes for each go process
#   genes_list <- list()
# 
#   # For each go process, retrieve the genes and store into genes_list
#   for (process in go_terms) {
#     go_term <- define_go_process[[process]]
#     if (is.null(go_term)) {
#       warning(paste("GO term not found for the process:", process))
#     } else {
#       genes <- get_genes_for_go_process(go_term)
#       genes_list[[process]] <- genes # Store list by process name
#     }
#   }
# 
#   print(genes_list)  # Debug print to check the expected structure
# 
#   genes_list  # Return the full list organized by process
# }



# Function to create the gene process mapping
create_gene_process_mapping <- function(genes_of_interest_for_heatmap_vector) {
  gene_process_mapping <- list()
  
  # Ensure correct indexing on the list
  for (process in names(genes_of_interest_for_heatmap_vector)) {
    genes <- genes_of_interest_for_heatmap_vector[[process]]  # ensure this accesses correct index
    for (gene in genes) {
      # Check symbol normalization for consistency
      gene <- toupper(gene)  # Normalize (if required) to match rownames format
      if (!is.null(gene_process_mapping[[gene]])) {
        gene_process_mapping[[gene]] <- c(gene_process_mapping[[gene]], process)
      } else {
        gene_process_mapping[[gene]] <- process
      }
    }
  }
  
  gene_process_mapping
}




# create_heatmap_combined_avg_log2fc.R
##############################################################################
# Prepare heatmap matrix should consider global variable dependency
prepare_heatmap_matrix <- function(data, genes_of_interest_for_heatmap_vector) {
  # Ensure both lists are uppercase
  genes_of_interest_for_heatmap_vector <- toupper(genes_of_interest_for_heatmap_vector)
  
  # Subset data
  subset <- data[rownames(data) %in% genes_of_interest_for_heatmap_vector, ]
  
  print("Subset after filtering:")
  print(subset)
  
  # Add process information
  subset$go_processes <- sapply(rownames(subset), function(gene) {
    gene <- toupper(gene)
    processes <- gene_process_mapping[[gene]]
    if (!is.null(processes)) {
      paste(processes, collapse = ",")
    } else {
      "None"
    }
  })
  
  as.matrix(subset[, 1:6])
}



# # Create heatmap function
# create_heatmap <- function(heatmap_matrix_numeric, go_association_vector) {
#   # Verify the input matrix
#   if (!is.matrix(heatmap_matrix_numeric)) {
#     stop("Input data must be a matrix.")
#   }
#   if (any(is.na(heatmap_matrix_numeric))) {
#     warning("NA values detected in the heatmap data.")
#   }
#   if (any(!is.finite(heatmap_matrix_numeric))) {
#     warning("Non-finite values detected in the heatmap data.")
#   }
#   
#   # If necessary, clean NA, NaN, Inf values
#   heatmap_matrix_numeric[is.na(heatmap_matrix_numeric)] <- 0
#   heatmap_matrix_numeric[!is.finite(heatmap_matrix_numeric)] <- 0
#   
#   # Define color function
#   col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
#   
#   # Create text annotation
#   row_annotation <- rowAnnotation(
#     text = anno_text(go_association_vector, gp = gpar(fontsize = 8)),
#     show_annotation_name = FALSE
#   )
#   
#   # Create the heatmap
#   heatmap <- Heatmap(
#     heatmap_matrix_numeric, 
#     name = "log2 fold change",
#     col = col_fun,
#     na_col = "grey",
#     cluster_rows = FALSE, 
#     show_row_names = TRUE,
#     row_names_side = "left",
#     row_names_gp = gpar(fontsize = 10),
#     row_gap = unit(3, "mm"),
#     cluster_columns = FALSE,  
#     show_column_names = TRUE,
#     column_names_rot = 45,
#     column_names_centered = TRUE,
#     column_title = "Differential Gene Expression AB vs SHAM",
#     show_column_dend = FALSE,
#     width = unit(10, "cm"),
#     heatmap_legend_param = list(
#       title = "avg_log2fc",
#       legend_direction = "vertical",
#       legend_height = unit(4, "cm")
#     ),
#     right_annotation = row_annotation,  # Use right_annotation, which is correct for adding text
#     gap = unit(1, "mm")
#   )
#   
#   # Create custom legend for NA
#   na_legend <- Legend(
#     labels = "NA",
#     title = "",
#     legend_gp = gpar(fill = "grey"),
#     direction = "vertical"
#   )
#   
#   # Draw the heatmap and stack legends
#   draw(heatmap, annotation_legend_list = list(na_legend), merge_legend = TRUE)
# }

# Create heatmap function with additional checks
create_heatmap <- function(heatmap_matrix_numeric, go_association_vector) {
  if (!is.matrix(heatmap_matrix_numeric)) {
    stop("Input data must be a matrix.")
  }
  
  heatmap_matrix_numeric[is.na(heatmap_matrix_numeric)] <- 0
  heatmap_matrix_numeric[!is.finite(heatmap_matrix_numeric)] <- 0
  
  if (nrow(heatmap_matrix_numeric) != length(go_association_vector)) {
    stop("Mismatch between number of matrix rows and association vector length.")
  }
  
  # Define color function
  col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  # Create text annotation with non-empty vector
  row_annotation <- rowAnnotation(
    text = anno_text(go_association_vector, gp = gpar(fontsize = 8)),
    show_annotation_name = FALSE
  )
  
  # Create the heatmap
  heatmap <- Heatmap(
    heatmap_matrix_numeric, 
    name = "log2 fold change",
    col = col_fun,
    na_col = "grey",
    cluster_rows = FALSE, 
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    row_gap = unit(3, "mm"),
    cluster_columns = FALSE,  
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_centered = TRUE,
    column_title = "Differential Gene Expression AB vs SHAM",
    show_column_dend = FALSE,
    width = unit(10, "cm"),
    heatmap_legend_param = list(
      title = "avg_log2fc",
      legend_direction = "vertical",
      legend_height = unit(4, "cm")
    ),
    right_annotation = row_annotation,  # Use right_annotation, which is correct for adding text
    gap = unit(1, "mm")
  )
  
  # Create custom legend for NA
  na_legend <- Legend(
    labels = "NA",
    title = "",
    legend_gp = gpar(fill = "grey"),
    direction = "vertical"
  )
  
  # Draw the heatmap and stack legends
  draw(heatmap, annotation_legend_list = list(na_legend), merge_legend = TRUE)
}

# ### Running script
# #######################################################################
# 
# #Load data
# mouse_vcm_all_genes_avg_log2fc <- readRDS("C:/Users/siljeew/snRNAseq/tmp/mouse_vcm_all_genes_avg_log2fc.rds")
# head(mouse_vcm_all_genes_avg_log2fc)
# 
# # # Use the row names from your dataset as gene names
# # gene_names <- rownames(mouse_vcm_all_genes_avg_log2fc)
# 
# # Define go processes of interest
# go_processes_of_interest_for_heatmap <- c('glycolysis')
# 
# # Retrieve genes of interest
# genes_of_interest_for_heatmap_vector <- get_genes_of_interest(go_processes_of_interest_for_heatmap, define_go_process, get_genes_for_go_process)
# View(genes_of_interest_for_heatmap_vector) # This contains glycolysis
# 
# # Process `genes_of_interest_for_heatmap_vector` to ensure it's flat
# genes_of_interest_for_heatmap_vector_unlisted <- unlist(genes_of_interest_for_heatmap_vector, use.names = FALSE)
# 
# # Create gene process mapping using the retrieved genes
# gene_process_mapping <- create_gene_process_mapping(genes_of_interest_for_heatmap_vector)
# 
# 
# # Prepare the heatmap matrix
# heatmap_matrix_numeric <- prepare_heatmap_matrix(mouse_vcm_all_genes_avg_log2fc, genes_of_interest_for_heatmap_vector_unlisted)
# nrow(heatmap_matrix_numeric) #71
# 
# # Check the column names and change to shorter and more readable names
# colnames(heatmap_matrix_numeric) <- c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")
# 
# # Create a vector for gene ontology mapping
# #######################################################################
# # Must be executed post `gene_process_mapping` calculation
# go_association_vector <- sapply(rownames(heatmap_matrix_numeric), function(gene) {
#   gene <- toupper(gene)  # Ensure consistent uppercase case
#   if (gene %in% names(gene_process_mapping)) {
#     processes <- gene_process_mapping[[gene]]
#     unique_processes <- unique(processes)
#     clean_process_name <- paste(unique_processes, collapse = ",")
#     clean_process_name <- gsub("[0-9]+", "", clean_process_name)
#     return(clean_process_name)
#   } else {
#     print(paste("No process found for gene:", gene))
#     return("None")
#   }
# }, USE.NAMES = FALSE)
# 
# # Check resulting vector
# print("Processed go_association_vector:")
# print(head(go_association_vector))
# 
# 
# ### Create and save heatmap
# #######################################################################
# heatmap_plot <- create_heatmap(heatmap_matrix_numeric, go_association_vector)
# 
# # Readable writing of names/legends/text?
# # Sort genes alphabetically after process?
