# Run create heatmap 
# Silje Wilhelmsen

# Source functions
########################################################################
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/load_libraries_run_create_heatmap_combined_avg_log2fc.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Analysis/define_go_terms_and_find_genes.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/gene_process_mapping.R")
source("C:/Users/siljeew/Master_project/snRNAseq/Visualization/create_heatmap_combined_avg_log2fc.R")

### Set parameters
#######################################################################
output_dir_heatmap <- "C:/Users/siljeew/Master_project/snRNAseq/Plots/Heatmap"

# Ensure the directory exists
if (!dir.exists(output_dir_heatmap)) {
  dir.create(output_dir_heatmap, recursive = TRUE)
}

#Load data
mouse_vcm_all_genes_avg_log2fc <- readRDS("C:/Users/siljeew/Master_project/snRNAseq/tmp/mouse_vcm_all_genes_avg_log2fc.rds")
head(mouse_vcm_all_genes_avg_log2fc)

# Define go processes of interest
go_processes_of_interest_for_heatmap <- c('transcription_coactivator_activity',
                                          'transcription_corepressor_activity',
                                          'dna_binding_transcription_activator_activity', 
                                          'dna_binding_transcription_repressor_activity',
                                          'ligand_modulated_transcription_factor_activity',
                                          'dna_binding_transcription_factor_activity_rna_polymeraseII_specific')


# Retrieve genes of interest
genes_of_interest_for_heatmap_vector <- get_genes_of_interest(go_processes_of_interest_for_heatmap, define_go_process, get_genes_for_go_process)
View(genes_of_interest_for_heatmap_vector)

# Process `genes_of_interest_for_heatmap_vector` to ensure it's flat
genes_of_interest_for_heatmap_vector_unlisted <- unlist(genes_of_interest_for_heatmap_vector, use.names = FALSE)

# Create gene process mapping using the retrieved genes
gene_process_mapping <- create_gene_process_mapping(genes_of_interest_for_heatmap_vector)


# Prepare the heatmap matrix
heatmap_matrix_numeric <- prepare_heatmap_matrix(mouse_vcm_all_genes_avg_log2fc, genes_of_interest_for_heatmap_vector_unlisted)
nrow(heatmap_matrix_numeric)

# Check the column names and change to shorter and more readable names
colnames(heatmap_matrix_numeric) <- c("6 Hours", "12 Hours", "1 Day", "3 Days", "1 Week", "3 Weeks")

# Create a vector for gene ontology mapping
#######################################################################
# Must be executed post `gene_process_mapping` calculation
go_association_vector <- sapply(rownames(heatmap_matrix_numeric), function(gene) {
  gene <- toupper(gene)  # Ensure consistent uppercase case
  if (gene %in% names(gene_process_mapping)) {
    processes <- gene_process_mapping[[gene]]
    unique_processes <- unique(processes)
    clean_process_name <- paste(unique_processes, collapse = ",")
    clean_process_name <- gsub("[0-9]+", "", clean_process_name)
    return(clean_process_name)
  } else {
    print(paste("No process found for gene:", gene))
    return("None")
  }
}, USE.NAMES = FALSE)

# Define color mapping for processes
process_colors <- c(
  'glycolysis' = "whitesmoke",
  'tca' = "wheat", 
  "glycolysis,tca" = "wheat4", 
  'pdh' = "slategray2",
  'glycolysis,pdh' = "thistle3",
  'fao' = "lightcoral",
  'etc' = "darkgreen",
  'chylomicron_clearance' = "pink", 
  'uncoupling_proteins' = "darkslategray3",
  'atp_synthesis' = "sienna2",
  'etc_succinate_to_ubiquinone' = "salmon3",
  'uncoupling_proteins' = "cornsilk2",
  'fatty_acid_transport' = "azure3",
  'b_oxidation' = "aliceblue",
  'monosaccharide_transmembrane_transport' = "papayawhip", 
  'protein_transmembrane_transport' = "lemonchiffon",
  'dna_binding_transcription_factor_activity' = "palegoldenrod",
  'translation_regulatory_activity' = "lightpink3",
  'ketone_metabolic_process' = "olivedrab",
  'hydroxybutyrate_dehydrogenase_activity' = "aquamarine4",
  'regulation_of_ketone_metabolic_process' = "orange2",
  'ketone_catabolism' = "rosybrown",
  'transcription_coactivator_activity' = "lemonchiffon",
  'transcription_corepressor_activity' = "whitesmoke",
  'dna_binding_transcription_activator_activity' = "papayawhip",
  'dna_binding_transcription_repressor_activity' = "slategray1",
  'ligand_modulated_transcription_factor_activity' = "aliceblue",
  'dna_binding_transcription_factor_activity_rna_polymeraseII_specific' = "cornsilk"
)



# Define a sorting order
sorting_order <- factor(go_association_vector, levels = c("transcription_coactivator_activity", 
                                                          "transcription_corepressor_activity", 
                                                          "dna_binding_transcription_activator_activity",
                                                          "dna_binding_transcription_repressor_activity",
                                                          "ligand_modulated_transcription_factor_activity",
                                                          "dna_binding_transcription_factor_activity_rna_polymeraseII_specific"))

# Order the matrix based on the sorting vector
ordered_indices <- order(sorting_order)

# Reorder heatmap matrix and process vector
heatmap_matrix_numeric <- heatmap_matrix_numeric[ordered_indices, ]
go_association_vector <- go_association_vector[ordered_indices]

# Verify sorting
print(head(rownames(heatmap_matrix_numeric)))
print(head(go_association_vector))

### Create and save heatmap
#######################################################################
heatmap_plot <- create_heatmap(heatmap_matrix_numeric, go_association_vector, process_colors)
any(is.na(heatmap_matrix_numeric))

save_heatmap_complex(heatmap_plot, output_dir_heatmap, go_processes_of_interest_for_heatmap)
