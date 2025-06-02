################################################################################
#####                   Multigroup heatmap for RNA-seq                     #####
################################################################################

# Create a single heatmap that includes genes enriched from several groups/comparisons

# Prior to running this script, you should have already run the appropriate 
# differential expression tests within Pluto 

# If the data has a number of equally interesting groups without a set of controls
# (e.g. a cancer dataset with several tumor locations), it may be desired to 
# compare each group to the rest of the dataset, such as bone vs. brain&lung&liver
# followed by brain vs. bone&lung&liver, etc.


################################################################################
#####                        User-defined parameters                       #####
################################################################################

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "xxx"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLXxxx"

# List the 'Plot ID' of differential expression results to include
plot_id <- c("xxxxxxx-xxx")
# Should be multiple Plot IDs for this analysis

# Which direction of gene changes to include?
fc_direction <- "positive" # Options: "positive", "negative", or "both".
# Assuming positive is most interesting as only genes with a positive fold-change
# in a comparison will be considered for inclusion
# If both is selected, the top `n_genes` in each direction will be included

# How many genes to consider for the heatmap **for each** comparison
n_genes <- 10

# Adjusted p-value cutoff to consider inclusion of a gene
significance_cutoff <- 0.05
# Once a gene meets this threshold, the selected genes are those with the largest 
# absolute fold-change based on `fc_direction`.

# Heatmap annotation?
annotation_column <- c(4,5) # The column number in the sample data to use for annotation
# Can be a single or multiple columns
# Set to NULL for no annotation
# Ideally, the annotation_column will relate to the groupings for the comparisons


################################################################################
#####                   Get genes to include in heatmap                    #####
################################################################################

# Load required libraries
library(pluto)
library(edgeR)
library(pheatmap)

# Log into Pluto
pluto_login(api_token)

# Loop each differential expression result and record included genes
merged_genes <- c()
for(result in plot_id){
  
  print(result)
  
  dge <- pluto_read_results(experiment_id, result)

  # Subset by `significance_cutoff`
  dge <- dge[dge$Adj_P_Value < significance_cutoff,]
  
  # Order by FC
  dge <- dge[order(dge$Log2_Fold_Change, decreasing = T),]
  
  # Top positive and negative genes
  dge_pos <- dge[dge$Log2_Fold_Change > 0,] # Ensure a positive FC
  genes_pos <- dge_pos$Gene_Symbol[1:n_genes]
  dge_neg <- dge[dge$Log2_Fold_Change < 0,] # Ensure a negative FC
  genes_neg <- dge_neg$Gene_Symbol[(nrow(dge_neg)-n_genes+1):nrow(dge_neg)]
  
  # Top genes
  if(fc_direction == "positive"){
    merged_genes <- c(merged_genes, genes_pos)
  }
  if(fc_direction == "negative"){
    merged_genes <- c(merged_genes, genes_neg)
  }
  if(fc_direction == "both"){
    merged_genes <- c(merged_genes, genes_pos, genes_neg)
  }
  
}
included_genes <- unique(merged_genes)
included_genes <- na.omit(included_genes)


################################################################################
#####             Extract normalized values of included genes              #####
################################################################################

# Read in experiment assay data
assay_data <- pluto_read_assay_data(experiment_id)

# Calculate counts per million (CPM) to get normalized data
if (any(duplicated(assay_data$gene_symbol))) {
  cat("Duplicate gene names detected in the first column\nDuplicated gene names will be made unique by appending a numeric suffix to the end (e.g. _1)")
  assay_data$gene_symbol <- make.unique(assay_data$gene_symbol, sep = "_")
}
rownames(assay_data) <- assay_data$gene_symbol
assay_data <- assay_data[,-1]
normalized_data <- data.frame(cpm(assay_data))
normalized_data <- log2(normalized_data + 1)

# Subset assay data to included genes
normalized_data <- normalized_data[rownames(normalized_data) %in% included_genes,]

# Z-score normalized_data
normalized_data <- data.frame(t(scale(t(normalized_data), center = TRUE, scale = TRUE)))


################################################################################
#####                          Generate heatmap                            #####
################################################################################

# Read in sample data
sample_data <- pluto_read_sample_data(experiment_id)

# Get annotation columns
if(!(is.null(annotation_column))){
  if(length(annotation_column) > 1){
    annotation_col <- data.frame(sample_data[,annotation_column])
    rownames(annotation_col) <- sample_data$sample_id
  } else {
    annotation_col <- data.frame(Group = sample_data[,annotation_column])
    colnames(annotation_col)[1] <- colnames(sample_data)[annotation_column]
    rownames(annotation_col) <- sample_data$sample_id
  }
}

# Generate heatmap
if(!(is.null(annotation_column))){
  
  ph <- pheatmap(normalized_data, scale = "none", 
           cluster_cols = TRUE, cluster_rows = TRUE,
           show_colnames = T,
           fontsize = 8,
           annotation_col = annotation_col)
  
} else{
    
    ph <- pheatmap(normalized_data, scale = "none", 
             cluster_cols = TRUE, cluster_rows = TRUE,
             show_colnames = T,
             fontsize = 8)
  
  }
ph

# Save plot
pdf("heatmap.pdf",
    width = 10, height = 6)
ph
dev.off()
