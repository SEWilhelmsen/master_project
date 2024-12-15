### Define GO terms and find related genes 
## Output: ??
# Silje Wilhelmsen

# Load libraries
library(AnnotationDbi)
library(org.Hs.eg.db)

# Define GO terms
define_go_process <- list(
  glycolysis = "GO:0006096",
  fao = "GO:0019395",
  tca = "GO:0006099"
)

get_genes_for_go_process <- function(define_go_process) {
  go_genes_for_process <- lapply(define_go_process, function(go_id) {
    entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, 
                                        keys = go_id, 
                                        keytype = "GOALL", 
                                        columns = "ENTREZID")[, "ENTREZID"]
    
    gene_symbols <- AnnotationDbi::select(org.Hs.eg.db, 
                                          keys = entrez_ids, 
                                          keytype = "ENTREZID", 
                                          columns = "SYMBOL")[, "SYMBOL"]
    gene_symbols
  })
  
  # Normalize gene symbols
  normalize_symbols <- function(genes) {
    toupper(genes)
  }
  
  normalized_go_genes <- lapply(go_genes_for_process, normalize_symbols)
  
  # Unlist go_genes and create a unique list of GO related genes
  all_go_genes <- unique(unlist(normalized_go_genes))  
  
  return(all_go_genes)
}