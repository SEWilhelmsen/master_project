### ANOVA
# Input: 
## Output: 
# Silje Wilhelmsen 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC6096346/ and https://statsandr.com/blog/anova-in-r/#preliminary-analyses 

# Preliminary analyses - visualize the data 

boxplot(expression ~ gene, 
        data = rna_data)


# Demultiplexing, Alignment, and Normalization
# Determining Intra- and Intergroup Sample Variability and Outliers
# - PCA, Pearsons correlation analysis, Spearmans correlation coefficient
# Count threshold?

# Identification of differential expressed genes (DEGs)
# - pairwise comparison between two groups and variance across groups