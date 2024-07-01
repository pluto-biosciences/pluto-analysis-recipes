################################################################################
#####    Marker expression analysis - single marker (cell scatter plot)    #####
################################################################################

# Script to create cell scatter plot visualizations for expression of a single
# marker gene for single cell RNA-seq (scRNA-seq) data in Pluto Bio.

# This example script illustrates how to create a cell scatter plot in UMAP
# space for a given marker gene. The cell scatter plot is equivalent to the
# Feature plot in Seurat (https://satijalab.org/seurat/), and can utilize t-SNE
# or PCA coordinates in place of UMAP coordinates.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

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
experiment_id <- "PLX141180"

# Define your target marker gene.
target <- "olig2"

# Define a file path for the analysis plot (cell scatter plot).
display_file_path <- paste0(experiment_id, "_", target, "_umap.html")

# Define a name for the analysis in Pluto.
analysis_name <- paste0("Marker expression analysis: ", target)

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Cell scatter plot using UMAP coordinates, colored by ",
    target, " expression."
)

# Colors for cell scatter plot color scale (low expression, high expression).
plot_colors <- c("#cdcfd1", "#AD05B1")

################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(tidyverse)
library(plotly)
library(htmlwidgets)

# Log into Pluto
pluto_login(api_token)

# Get Seurat object
exp_obj <- pluto_read_seurat_object(
    experiment_id = experiment_id,
    seurat_type = "final"
)
so <- exp_obj$obj

# Plot target expression in UMAP space
# Note: to plot target expression in t-SNE space, replace instances of 'umap'/
# 'UMAP' with 'tsne'/'tSNE'. To plot target expression in PCA space, replace
# instances of 'umap'/'UMAP' with 'pca'/'PC'.

# Get plot data (target expression and coordinates)
plot_data <- bind_cols(
    as.data.frame(so@reductions$umap@cell.embeddings),
    as.data.frame(so@assays$RNA@data[`target`, ])
)
colnames(plot_data)[ncol(plot_data)] <- target

# Create cell scatter plot
cell_scatter_plot <- plot_ly(plot_data) %>%
    add_trace(
        x = plot_data$UMAP_1,
        y = plot_data$UMAP_2,
        type = "scattergl",
        mode = "markers",
        color = plot_data[[target]],
        colors = plot_colors,
        marker = list(size = 3),
        legendgroup = plot_data[[target]],
        showlegend = TRUE
    ) %>%
    colorbar(title = "Expr.") %>%
    layout(
        xaxis = list(
            title = list(text = "UMAP 1"),
            showgrid = FALSE,
            zeroline = FALSE,
            showline = TRUE
        ),
        yaxis = list(
            title = list(text = "UMAP 2"),
            showgrid = FALSE,
            zeroline = FALSE,
            showline = TRUE
        )
    ) %>%
    config(
        displayModeBar = FALSE
    ) %>%
    style(hoverinfo = "none")

# Save cell scatter plot
saveWidget(
    cell_scatter_plot,
    file = display_file_path,
    title = paste0("Cell scatter plot, colored by ", target, " expression"),
    selfcontained = TRUE
)

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)
