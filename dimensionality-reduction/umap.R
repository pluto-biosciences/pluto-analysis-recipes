################################################################################
#####              Uniform manifold approximation & projection             #####
################################################################################

# Script to map samples into a low-dimensional space with uniform manifold
# approximation & projection (UMAP) to examine patterns in a Pluto Bio
# experiment.

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
# Note: at least 2 samples are required to run UMAP.
# sample_ids <- c("R6_rep1", "R6_rep2", "R6_WT_H3_3_rep1", "R6_WT_H3_3_rep2")
sample_ids <- NULL # This will include all samples

# Define the number of nearest neighbors to use (maximum of 12).
n_neighbors <- 3

# Define the minimum distance between points.
min_distance <- 0.5

# Perform CPM-normalization on the assay data?
# Set this to TRUE if working with raw counts data!
cpm_normalize <- TRUE

# Define variables to group samples by.
# Note: these must match column names in your sample data!
group_variables <- c("cell_line_modification")

# Define a name for the analysis in Pluto.
analysis_name <- "Uniform manifold approximation & projection"

# Define a file path for the analysis plot (scatter plot).
display_file_path <- paste0(experiment_id, "_umap_scatter_plot.html")

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(experiment_id, "_umap_results.csv")

# Include methods to describe your analysis.
if (is.null(sample_ids)) {
    plot_methods <- paste0(
        "Uniform manifold approximation & projection (UMAP) was performed by ",
        "applying the umap() R funciton to counts per million (CPM)-",
        "normalized values for all targets in the experiment and all samples. ",
        "The minimum point distance used was ", min_distance, " and the ",
        "number of neighbors was set to ", n_neighbors, "."
    )
} else {
    plot_methods <- paste0(
        "Uniform manifold approximation & projection (UMAP) was performed by ",
        "applying the umap() R funciton to counts per million (CPM)-",
        "normalized values for all targets in the experiment and the ",
        "following sample IDs: ",
        paste(sample_ids, collapse = ", "), ". The minimum point distance ",
        "used was ", min_distance, " and the number of neighbors was set to ",
        n_neighbors, "."
    )
}

################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(edgeR)
library(umap)
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
#----                               Run UMAP                               ----#
#------------------------------------------------------------------------------#
# Create custom UMAP analysis config to store desired params
custom_config <- umap.defaults
custom_config$n_neighbors <- n_neighbors
custom_config$min_dist <- min_distance
custom_config$n_components <- min(c(5, ncol(assay_data)))
custom_config$random_state <- 7
custom_config$transform_state <- 42
# print(custom_config)

# Transpose assay data
t_assay_data <- as.data.frame(t(assay_data))

# Only include genes with non-zero variance
to_include <- unlist(lapply(t_assay_data, function(l) {
    v <- var(l)
    return(!is.na(v) & v > 0)
}))
t_assay_data <- t_assay_data[, to_include]

# Run UMAP with custom config
assay_data_umap <- umap(t_assay_data, config = custom_config)
umap_table <- as.data.frame(assay_data_umap$layout)
names(umap_table) <- paste0("UMAP", seq_len(ncol(umap_table)))
umap_table <- cbind(sample_id = row.names(umap_table), umap_table)

# Save results table
write.csv(umap_table, results_file_path, row.names = FALSE)

#------------------------------------------------------------------------------#
#----                         Create a scatter plot                        ----#
#------------------------------------------------------------------------------#
# Create a scatter plot, Pluto style!
plot_data <- merge(sample_data, umap_table, by = "sample_id", all.x = TRUE)

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
    x = ~UMAP1,
    y = ~UMAP2,
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
            title = "UMAP1",
            showgrid = FALSE,
            zeroline = TRUE,
            zerolinecolor = "lightgray",
            showline = TRUE,
            linecolor = "black"
        ),
        yaxis = list(
            title = "UMAP2",
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
    title = "Uniform manifold approximation & projection",
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
