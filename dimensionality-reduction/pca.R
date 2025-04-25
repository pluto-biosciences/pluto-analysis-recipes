################################################################################
#####                     Principal components analysis                    #####
################################################################################

# Script to map samples into a low-dimensional space with principal components
# analysis (PCA) to examine patterns in a Pluto Bio experiment.

# Note that this script is designed for bulk data (i.e., bulk RNA-seq, ATAC-seq,
# ChIP-seq), not single cell / nuclei data.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else,
# with the exception of adjusting the plot colors in plotting portion of the
# Main script.

# For more advanced plot customization, scroll down to the Main script and edit
# parameters for the  plot_ly function. For more information, see the usage
# docs: https://plotly.com/r/line-and-scatter/

################################################################################
#####                        User-defined parameters                       #####
################################################################################

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX271798"

# List the sample IDs you want to use for the analysis. Set this field to NULL
# to use all samples in the experiment.
# Note: at least 2 samples are required to run PCA.
# sample_ids <- c("R6_rep1", "R6_rep2", "R6_WT_H3_3_rep1", "R6_WT_H3_3_rep2")
sample_ids <- NULL # This will include all samples

# Shift data to be zero centered?
center_data <- TRUE

# Scale data to have unit variance?
scale_data <- TRUE

# Perform CPM-normalization on the assay data?
# Set this to TRUE if working with raw counts data!
cpm_normalize <- TRUE

# Define variables to group samples by.
# Note: these must match column names in your sample data!
group_variables <- c("cell_line_modification")

# Define whether to use variable features for PCA (optional).
use_var_features <- FALSE # Set to TRUE to use variable features for PCA
num_var_features <- 1000 # Number of variable features to use for PCA

# Define whether to use batch adjustment (optional).
batch_adjustment <- FALSE # Set to TRUE to perform batch adjustment with limma
batch_adjustment_var <- "batch" # Column name in sample data for batch adjustment

# Define a name for the analysis in Pluto.
analysis_name <- "Principal components analysis"

# Define a file path for the analysis plot (scatter plot).
display_file_path <- paste0(experiment_id, "_pca_scatter_plot.html")

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(experiment_id, "_pca_results.csv")

#------------------------------------------------------------------------------#
#----                                Methods                               ----#
#------------------------------------------------------------------------------#
# Include methods to describe your analysis.
plot_methods <- paste0(
    "Principal components analysis (PCA) was calculated by applying the ",
    "prcomp() R function to "
)
if (cpm_normalize) {
    plot_methods <- paste0(
        plot_methods,
        "counts per million (CPM)-normalized values for all targets in the ",
        "experiment and "
    )
} else {
    plot_methods <- paste0(
        plot_methods,
        "raw counts values for all targets in the experiment and "
    )
}
if (is.null(sample_ids)) {
    plot_methods <- paste0(
        plot_methods,
        "all samples."
    )
} else {
    plot_methods <- paste0(
        plot_methods,
        "the following sample IDs: ",
        paste(sample_ids, collapse = ", "), "."
    )
}
if (batch_adjustment) {
    plot_methods <- paste0(
        plot_methods,
        " Normalized data was batch-corrected for ", batch_adjustment_var,
        " using the removeBatchEffect() function in limma."
    )
}
if (center_data) {
    plot_methods <- paste0(
        plot_methods,
        " The data was shifted to be zero centered."
    )
}
if (scale_data) {
    plot_methods <- paste0(
        plot_methods,
        " The data was scaled to have unit variance before PCA was computed."
    )
}
if (use_var_features) {
    plot_methods <- paste0(
        plot_methods,
        " The top ", num_var_features, " features with the highest variance ",
        "were used across "
    )
    if (is.null(sample_ids)) {
        plot_methods <- paste0(
            plot_methods,
            "all samples."
        )
    } else {
        plot_methods <- paste0(
            plot_methods,
            "the included samples."
        )
    }
}

################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(edgeR)
library(limma)
library(tidyverse)
library(plotly)
library(htmlwidgets)

#------------------------------------------------------------------------------#
#----                             Prepare data                             ----#
#------------------------------------------------------------------------------#
# Log into Pluto
pluto_login(api_token)

# Read in experiment assay data
assay_data <- pluto_read_assay_data(experiment_id)

# Use the first column (gene symbol or peak ID) for the features
row.names(assay_data) <- make.names(assay_data[, 1], unique = TRUE)
assay_data[, 1] <- NULL

# Read in experiment sample data
sample_data <- pluto_read_sample_data(experiment_id)

# Subset assay_data and sample_data to selected sample IDs, or include all
# samples
if (!is.null(sample_ids)) {
    assay_data <- assay_data[, names(assay_data) %in% sample_ids]
    sample_data <- sample_data %>% filter(sample_id %in% sample_ids)
} else { # include all samples
    # This will drop any extra, non-sample columns (gene_symbol, probe_id, etc.)
    assay_data <- assay_data[, colnames(assay_data) %in% sample_data$sample_id]
}

# Cast all columns to numeric values
if (!all(unlist(lapply(assay_data, is.numeric)))) {
    nonnumeric_columns <- names(assay_data)[which(
        !unlist(lapply(assay_data, is.numeric))
    )]
    assay_data <- as.data.frame(lapply(assay_data, as.numeric))
}
message(
    "Assay data contains ", nrow(assay_data), " features and ",
    ncol(assay_data), " samples total"
)

# Perform CPM-normalization
if (cpm_normalize) {
    assay_data <- as.data.frame(edgeR::cpm(assay_data, log = FALSE))
    assay_data <- log2(assay_data + 1)
}

#------------------------------------------------------------------------------#
#----                           Batch correction                           ----#
#------------------------------------------------------------------------------#
# Perform limma batch correction if batch_adjustment is TRUE.
# Note: batch adjustment with limma must be performed on normalized data.
if (cpm_normalize && batch_adjustment) {
    # Check that batch column has more than one value
    if (length(unique(sample_data[[batch_adjustment_var]])) < 2) {
        message("Batch adjustment column has fewer than 2 unique values.")
        message("Skipping batch adjustment...")
    }
    # Batch adjustment with limma
    assay_data <- data.frame(
        removeBatchEffect(
            assay_data,
            batch = sample_data[[batch_adjustment_var]]
        )
    )
} else if (!cpm_normalize && batch_adjustment) {
    message("Batch adjustment with limma requires CPM-normalized data.")
    message("Skipping batch adjustment...")
}

#------------------------------------------------------------------------------#
#----                   Subset to most variable features                   ----#
#------------------------------------------------------------------------------#
if (use_var_features) {
    if (num_var_features > nrow(assay_data)) {
        num_var_features <- nrow(assay_data)
        message("num_var_features > number of rows in assay_data.")
        message("Defaulting to ", num_var_features, " features.")
    }
    message("Subsetting to top ", num_var_features, " variable features...")
    sample_variance <- apply(assay_data, 1, var)
    sample_variance_ordered <- sort(sample_variance, decreasing = TRUE)
    assay_data <- assay_data[
        rownames(assay_data) %in% names(
            sample_variance_ordered
        )[1:num_var_features],
    ]
}

#------------------------------------------------------------------------------#
#----                                Run PCA                               ----#
#------------------------------------------------------------------------------#
# Transpose assay data
t_assay_data <- as.data.frame(t(assay_data))

# Only include genes with non-zero variance
to_include <- unlist(lapply(t_assay_data, function(l) {
    v <- var(l)
    return(!is.na(v) & v > 0)
}))
t_assay_data <- t_assay_data[, to_include]

# Run PCA
pca_results <- prcomp(
    x = t_assay_data, center = center_data, scale. = scale_data
)
pca_table <- as.data.frame(pca_results$x)

# Limit to 5 PC columns max
if (ncol(pca_table) > 5) {
    pca_table <- pca_table[, 1:5]
}

# Add sample ID column
pca_table <- cbind(sample_id = row.names(pca_table), pca_table)

# Save results table
write.csv(pca_table, results_file_path, row.names = FALSE)

#------------------------------------------------------------------------------#
#----                         Create a scatter plot                        ----#
#------------------------------------------------------------------------------#
# Create a scatter plot, Pluto style!
plot_data <- merge(sample_data, pca_table, by = "sample_id", all.x = TRUE)

# Add group info
if (length(group_variables) == 1) {
    plot_data$color_by <- plot_data[[group_variables]]
} else {
    plot_data$color_by <- apply(
        plot_data[, group_variables], 1, function(row) {
            paste(
                row,
                collapse = " "
            )
        }
    )
}

# Add custom colors for each group
groups <- unique(plot_data$color_by)
# print(groups)
custom_colors <- c(
    "#312E81", "#FCA5A5", "#84CC16", "#F59E0B", "#0E7490", "#F5D0FE"
)
color_mapping <- setNames(custom_colors, groups)
plot_data$color <- color_mapping[plot_data$color_by]

# Create the scatter plot
scatter_plot <- plot_ly(
    data = plot_data,
    x = ~PC1,
    y = ~PC2,
    type = "scatter",
    mode = "markers",
    text = ~sample_id, # Hover text with sample ID
    hoverinfo = "text", # Show only sample ID on hover
    hovertemplate = paste(
        "%{text}",
        "<br>x: %{x}",
        "<br>y: %{y}",
        "<extra></extra>"
    ),
    color = ~color_by,
    colors = ~color,
    marker = list(
        size = 11, # Point size
        opacity = 1 # Point opacity
    )
) %>%
    layout(
        xaxis = list(
            title = "PC1",
            showgrid = FALSE,
            zeroline = TRUE,
            zerolinecolor = "lightgray",
            showline = TRUE,
            linecolor = "black"
        ),
        yaxis = list(
            title = "PC2",
            showgrid = FALSE,
            zeroline = TRUE,
            zerolinecolor = "lightgray",
            showline = TRUE,
            linecolor = "black"
        ),
        showlegend = TRUE
    )

# Save scatter plot
saveWidget(
    scatter_plot,
    file = display_file_path,
    title = "Principal components analysis",
    selfcontained = TRUE
)

#------------------------------------------------------------------------------#
#----                         Push results to Pluto                        ----#
#------------------------------------------------------------------------------#
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    results_file_path = results_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)
