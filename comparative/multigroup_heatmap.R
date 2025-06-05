################################################################################
#####               Multigroup heatmap for bulk RNA-seq data               #####
################################################################################

# Create a single heatmap that includes genes enriched from several groups or
# comparisons.

# Prior to running this script, you should have already run the appropriate
# differential expression tests within your Pluto Bio experiment.

# Note that this script is designed for bulk RNA-seq data.

# If the data has a number of equally interesting groups without a set of
# controls (e.g. a cancer dataset with several tumor locations), it may be
# desired to compare each group to the rest of the dataset, such as bone vs
# brain&lung&liver followed by brain vs bone&lung&liver, etc.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

# For more advanced heatmap customization, scroll down to the Main script and
# edit parameters for the ComplexHeatmap Heatmap function. For more information,
# see the usage docs: https://jokergoo.github.io/ComplexHeatmap-reference/book/

################################################################################
#####                        User-defined parameters                       #####
################################################################################

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX162695"

# List the plot IDs of differential expression results to include
plot_ids <- c(
  "821ae6bb-fde4-4ca2-98c1-f9382b23d7f5",
  "fd34fa5a-ba3f-4036-a621-4ed822b858a4",
  "d5e64986-b283-497b-8d60-fab14c87ed10"
)

# Which direction of gene changes to include?
# i.e., upregulated (positive), downregulated (negative), or both.
# If "both" is selected, the top `n_genes` in each direction are included.
fc_direction <- "positive" # Options: "positive", "negative", or "both".

# How many genes for each comparison should be included in the heatmap?
n_genes <- 20

# What adjusted p-value cutoff should be used to select genes?
# This is the threshold for significance in the differential expression results.
# Note that selected genes are those that meet this threshold and have the
# largest absolute fold-change based on `fc_direction`.
significance_cutoff <- 0.01

# Define the column number in the sample data to use for annotation.
# This can be a single column or multiple columns. Set to NULL for no
# annotation. Ideally, the annotation column will relate to the groupings used
# for the comparisons.
annotation_column <- c(2, 5)

# Define a file path for the analysis plot (heatmap) (PDF only).
display_file_path <- paste0(experiment_id, "_multigroup_heatmap.pdf")

# Define a name for the analysis in Pluto.
analysis_name <- "Multigroup heatmap of top genes"

# Include methods to describe the analysis.
plot_methods <- paste0(
  "Differential expression results were filtered by adjusted p-value ",
  "and ranked by log2 fold change. For each comparison, the top upregulated ",
  "and/or downregulated genes were selected. Gene expression values were ",
  "CPM-normalized, log2-transformed, and then z-scored by row to standardize ",
  "expression across samples. The resulting matrix was visualized as a ",
  "heatmap using the pheatmap R package, with rows representing genes and ",
  "columns representing samples, optionally annotated by sample metadata."
)

################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(edgeR)
library(ComplexHeatmap)
library(circlize)

# Log into Pluto
pluto_login(api_token)

#------------------------------------------------------------------------------#
#----                   Get top genes for each comparison                  ----#
#------------------------------------------------------------------------------#
top_genes <- c()
for (result in plot_ids) {
  dge <- pluto_read_results(experiment_id, result)

  # Subset by `significance_cutoff`
  dge <- dge[dge$Adj_P_Value < significance_cutoff, ]

  # Order by FC
  dge <- dge[order(dge$Log2_Fold_Change, decreasing = TRUE), ]

  # Top positive and negative genes
  dge_pos <- dge[dge$Log2_Fold_Change > 0, ] # Ensure a positive FC
  genes_pos <- dge_pos$Gene_Symbol[1:n_genes]
  dge_neg <- dge[dge$Log2_Fold_Change < 0, ] # Ensure a negative FC
  genes_neg <- dge_neg$Gene_Symbol[(nrow(dge_neg) - n_genes + 1):nrow(dge_neg)]

  # Top genes
  if (fc_direction == "positive") {
    top_genes <- c(top_genes, genes_pos)
  }
  if (fc_direction == "negative") {
    top_genes <- c(top_genes, genes_neg)
  }
  if (fc_direction == "both") {
    top_genes <- c(top_genes, genes_pos, genes_neg)
  }
}
top_genes <- unique(top_genes)
top_genes <- na.omit(top_genes)

#------------------------------------------------------------------------------#
#----                Extract normalized values of top genes                ----#
#------------------------------------------------------------------------------#
# Read in experiment assay data
assay_data <- pluto_read_assay_data(experiment_id)

# Calculate counts per million (CPM) to get normalized data
if (any(duplicated(assay_data$gene_symbol))) {
  cat(
    "Duplicate gene names detected in the first column\n",
    "Duplicated gene names will be made unique by appending a numeric suffix ",
    "to the end (e.g. _1)"
  )
  assay_data$gene_symbol <- make.unique(assay_data$gene_symbol, sep = "_")
}
rownames(assay_data) <- assay_data$gene_symbol
assay_data <- assay_data[, -1]
normalized_data <- data.frame(cpm(assay_data))
normalized_data <- log2(normalized_data + 1)

# Subset assay data to top genes
normalized_data <- normalized_data[rownames(normalized_data) %in% top_genes, ]

# Z-score normalized_data
normalized_data <- data.frame(
  t(scale(t(normalized_data), center = TRUE, scale = TRUE))
)
normalized_data_mat <- as.matrix(normalized_data)

#------------------------------------------------------------------------------#
#----                    Prepare multigroup heatmap data                   ----#
#------------------------------------------------------------------------------#
# Get heatmap annotations
# ***** BE SURE TO ADJUST THE ANNOTATION COLORS FOR YOUR ANALYSIS!! *****
if (!is.null(annotation_column)) {
  sample_data <- pluto_read_sample_data(experiment_id)
  matched_samples <- colnames(normalized_data_mat)
  sample_data <- sample_data[sample_data$sample_id %in% matched_samples, ]
  sample_data <- sample_data[match(matched_samples, sample_data$sample_id), ]

  # Create annotation_col
  if (length(annotation_column) > 1) {
    annotation_col <- data.frame(sample_data[, annotation_column, drop = FALSE])
  } else {
    annotation_col <- data.frame(Group = sample_data[, annotation_column])
    colnames(annotation_col)[1] <- colnames(sample_data)[annotation_column]
  }

  # Set rownames for annotation to match column names of normalized data
  rownames(annotation_col) <- sample_data$sample_id

  # Define annotation colors for each annotation column
  # Note: adjust this manually for your specific data!!!
  group_color_pal <- c("#FCA5A5", "#0E7490", "#65A30D", "#F59E0B")
  group_levels <- unique(annotation_col$group)
  group_colors <- setNames(
    group_color_pal, group_levels
  )

  tissue_color_pal <- c("#312e81", "#F5D0FE")
  tissue_levels <- unique(annotation_col$tissue)
  tissue_colors <- setNames(
    tissue_color_pal, tissue_levels
  )

  column_ha <- HeatmapAnnotation(
    df = annotation_col,
    col = list(
      group = group_colors,
      tissue = tissue_colors
    )
  )
}

# Create color function for heatmap values
col_fun <- colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

#------------------------------------------------------------------------------#
#----                      Generate multigroup heatmap                     ----#
#------------------------------------------------------------------------------#
# Create heatmap
# ***** CUSTOMIZE THE HEATMAP() FUNCTION OPTIONS FOR YOUR ANALYSIS!! *****
if (!(is.null(annotation_column))) {
  hm <- Heatmap(
    normalized_data_mat,
    name = "Z-score",
    col = col_fun,
    cluster_rows = TRUE, # add row dendrogram?
    cluster_columns = FALSE, # add column dendrogram?
    show_column_names = TRUE, # show column (sample) names?
    show_row_names = TRUE, # show row (gene) names?
    row_names_gp = gpar(fontsize = 8), # row name fontsize
    column_names_gp = gpar(fontsize = 8), # column name fontsize
    # column_split = annotation_col$group, # split columns by annotation column
    top_annotation = column_ha
  )
} else {
  hm <- Heatmap(
    normalized_data_mat,
    name = "Z-score",
    col = col_fun,
    cluster_rows = TRUE, # add row dendrogram?
    cluster_columns = TRUE, # add column dendrogram?
    show_column_names = TRUE, # show column (sample) names?
    row_names_gp = gpar(fontsize = 8), # row name fontsize
    column_names_gp = gpar(fontsize = 8), # column name fontsize
    show_row_names = TRUE, # show row (gene) names?
  )
}
hm

# Save the heatmap to a PDF file
pdf(display_file_path, width = 7, height = 7)
hm
dev.off()

#------------------------------------------------------------------------------#
#----                         Push results to Pluto                        ----#
#------------------------------------------------------------------------------#
pluto_add_experiment_plot(
  experiment_id = experiment_id,
  display_file_path = display_file_path,
  analysis_name = analysis_name,
  plot_methods = plot_methods
)
