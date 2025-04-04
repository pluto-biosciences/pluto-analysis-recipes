################################################################################
#####                          Clustering analysis                         #####
################################################################################

# Script to perform clustering analysis on genome-wide feature profiles and
# create a heatmap visualization using the ComplexHeatmap R package.

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
experiment_id <- "PLX203595"

# Provide a CSV file containing your features (genes) of interest.
# The file should contain a single column with the header "gene_symbol", and an
# optional "category" column if you want to group your heatmap rows.
# Example file structure:
#     gene_symbol,category
#     POU5F1,Pluripotency
#     NANOG,Pluripotency
#     SOX2,Neural ectoderm
#     SOX21,Neural ectoderm
gene_list <- "nociceptor_markers.csv"

# Optional: cluster rows (genes) by their category in the gene_list CSV file.
# To ignore this parameter, set it to FALSE: group_by_category <- FALSE
group_by_category <- TRUE

# Optional: define an order for the gene categories.
# This is useful if you want to set the order of the heatmap rows.
# To ignore this parameter, set it to NULL: category_levels <- NULL
category_levels <- c(
    "Pluripotency", "Neural ectoderm", "Neural crest",
    "Nociceptor", "Peptidergic", "Nonpeptidergic", "Neurofilament"
)

# List the sample IDs you want to use for the analysis.
# These must match the sample IDs in the experiment sample data!
# Example: samples <- c("Nociceptors_D4_rep1", "Nociceptors_D8_rep1")
# To use all samples, set this parameter to NULL.
samples <- NULL

# Optional: define a variable to group samples by and annotate the heatmap.
# The grouping column must match a variable in the experiment sample data!
# To ignore this parameter, set it to NULL: grouping_col <- NULL
grouping_col <- "age"

# Optional: define an order for the grouping variable values.
# This is useful if you want to set the order of the heatmap columns.
# To ignore this parameter, set it to NULL: grouping_levels <- NULL
grouping_levels <- c("0 days", "4 days", "8 days", "12 days", "28 days")

# Optional: define a color scheme for the column annotations. You need one
# color for each value in the grouping_levels parameter listed above.
# To ignore this parameter, set it to NULL: annot_colors <- NULL
group_colors <- c("#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51")

# Define a file path for the analysis plot (heatmap) (PDF only).
display_file_path <- paste0(experiment_id, "_clustering_heatmap.pdf")

# Include methods to describe the analysis.
plot_methods <- paste0(
    "CPM-normalized gene expression data was z-scored by row (i.e., each ",
    "gene's expression was normalized across samples) to standardize the ",
    "values and make them comparable across genes. The resulting data was ",
    "then visualized as a heatmap using the ComplexHeatmap R package, where ",
    "rows represent genes and columns represent samples."
)

################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(edgeR)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

#------------------------------------------------------------------------------#
#----                             Prepare data                             ----#
#------------------------------------------------------------------------------#
# Log into Pluto
pluto_login(api_token)

# Read in experiment assay data
assay_data <- pluto_read_assay_data(experiment_id)

# Read in experiment sample data
sample_data <- pluto_read_sample_data(experiment_id)

# Filter assay and sample data if samples is provided
if (!is.null(samples)) {
    assay_data <- assay_data[, colnames(assay_data) %in% samples]
    sample_data <- subset(sample_data, sample_id %in% samples)
}

# Perform CPM-normalization on assay data
rownames(assay_data) <- assay_data$gene_symbol
assay_data <- assay_data[, -1]
assay_data_cpm <- data.frame(cpm(assay_data))
assay_data_cpm <- log2(assay_data_cpm + 1)

# Prepare the column (sample) annotations if grouping_col is provided
if (!is.null(grouping_col)) {
    annot_col <- as.data.frame(sample_data[[grouping_col]])
    row.names(annot_col) <- sample_data$sample_id
    colnames(annot_col)[1] <- grouping_col
    annot_col$sample_id <- sample_data$sample_id
    if (!is.null(grouping_levels)) { # Order groups by grouping_levels
        annot_col[[grouping_col]] <- factor(
            annot_col[[grouping_col]],
            levels = grouping_levels
        )
        annot_col <- annot_col[order(annot_col[[grouping_col]]), ]
        sample_order <- annot_col$sample_id
        annot_col$sample_id <- NULL
        assay_data_cpm <- assay_data_cpm[, sample_order]
    }

    # Define column (sample) annotation colors
    if (!is.null(group_colors)) {
        annot_col_colors <- list(
            setNames(group_colors, grouping_levels)
        )
        names(annot_col_colors) <- grouping_col
    } else {
        color_pal <- brewer.pal(n = length(grouping_levels), name = "Spectral")
        annot_col_colors <- list(
            setNames(color_pal, grouping_levels)
        )
        names(annot_col_colors) <- grouping_col
    }
}

# Subset assay data to features (genes) of interest
assay_data_matrix <- as.matrix(assay_data_cpm)
genes <- read.csv(gene_list)
matched_indices <- match(genes$gene_symbol, rownames(assay_data_matrix))
genes_in_matrix <- genes$gene_symbol[!is.na(matched_indices)]
plt_matrix <- assay_data_matrix[genes_in_matrix, ]

# Perform row Z-score normalization
cal_z_score <- function(x) {
    (x - mean(x)) / sd(x)
}
plt_matrix_norm <- t(apply(plt_matrix, 1, cal_z_score))
plt_matrix_norm <- plt_matrix_norm[complete.cases(plt_matrix_norm), ]

# Prepare the row (gene) annotations if group_by_category is set to TRUE
if (group_by_category) {
    gene_col <- genes
    rownames(gene_col) <- genes$gene_symbol
    gene_col$gene_symbol <- NULL
    if (!is.null(category_levels)) { # Order categories by category_levels
        gene_col$category <- factor(
            gene_col$category,
            levels = category_levels
        )
    }
}

#------------------------------------------------------------------------------#
#----                           Create a heatmap                           ----#
#------------------------------------------------------------------------------#
# Define the heatmap color scale
heatmap_colors <- colorRamp2(
    seq(-3, 3, length = 20),
    colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(20)
)

# Set cluster_rows based on group_by_category
if (group_by_category) {
    cluster_rows <- FALSE
} else {
    cluster_rows <- TRUE
}

# Set cluster_columns and show_column_names based on grouping_col
if (is.null(grouping_col)) {
    cluster_columns <- TRUE
    show_column_names <- TRUE
} else {
    cluster_columns <- FALSE
    show_column_names <- FALSE
}

# Conditionally create top_annotation only if grouping_col is not NULL
if (!is.null(grouping_col)) {
    top_annotation <- HeatmapAnnotation(
        df = annot_col, # Data frame with annotation info
        col = annot_col_colors, # Annotation colors
        annotation_name_side = "left"
    )
} else {
    top_annotation <- NULL
}

# Conditionally set row_split, row_title_rot, and row_title_gp if
# group_by_category is TRUE
if (group_by_category) {
    row_split <- gene_col
} else {
    row_split <- NULL
}

# Create the heatmap using ComplexHeatmap
heatmap <- Heatmap(
    plt_matrix_norm,
    name = "Z-score", # Heatmap legend title
    col = heatmap_colors, # Heatmap color palette
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_column_names = show_column_names,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8), # Row name font size
    column_names_gp = gpar(fontsize = 8), # Column name font size
    top_annotation = top_annotation,
    row_split = row_split,
    row_title_rot = 0,
    row_title_gp = gpar(fontsize = 10)
)

# Save the heatmap to a PDF file
pdf(display_file_path, width = 7, height = 5)
draw(heatmap)
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
