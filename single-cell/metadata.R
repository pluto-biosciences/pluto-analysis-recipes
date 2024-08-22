################################################################################
#####                           Metadata analysis                          #####
################################################################################

# Script to create cell scatter plot visualization, grouping cells by a metadata 
# variable(s), for single cell RNA-seq (scRNA-seq) data in Pluto Bio.

################################################################################
#####                           Cell scatter plot                          #####
################################################################################

# This example script illustrates how to create a cell scatter plot in UMAP
# space for a given metadata variable(s). The cell scatter plot is equivalent to 
# the DimPlot in Seurat (https://satijalab.org/seurat/), and can utilize t-SNE
# or PCA coordinates in place of UMAP coordinates.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

# For more advanced plot customization, scroll down to the Main script (Cell
# scatter plot) and edit parameters for the  plot_ly function. For more
# information, see the usage docs: https://plotly.com/r/line-and-scatter/

#------------------------------------------------------------------------------#
#----                        User-defined parameters                       ----#
#------------------------------------------------------------------------------#

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX287028"

# Define your grouping variable(s).
# Remember, you can access available grouping variables by running:
# head(so@meta.data)
cell_grouping <- c("default_clusters_res_0_2")

# Define a file path for the analysis plot (cell scatter plot).
display_file_path <- paste0(
    experiment_id, "_umap_grouped_by_", 
    paste(cell_grouping, collapse = "_"), 
    ".html"
)

# Define a name for the analysis in Pluto.
analysis_name <- paste0(
    "Metadata analysis: ", 
    paste(cell_grouping, collapse = ", ")
)

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Cell scatter plot using UMAP coordinates, colored by ",
    paste(cell_grouping, collapse = ", ")
)

# Colors for cell groups.
# If you are grouping cells by one of your clustering resolutions, you can
# access the same color palette you are using in the Pluto app by loading in
# your custom color palettes and selecting the palette that corresponds to your
# resolution of interest after loading in exp_obj (see Main script).
custom_palettes <- exp_obj$colors
plot_colors <- custom_palettes$default_clusters_res_0_2
# For more info, see:
# https://pluto-biosciences.github.io/pluto-sdk-R/articles/scrnaseq_recipes.html

# You can also group your cells by multiple metadata variables. For example:
# cell_grouping <- c("sex", "genotype")
#
# Note: the length of plot_colors must match the number of combinatorial
# categories contained in cell_grouping. In the above example, the metadata 
# variable "sex" contains "Female" and "Male", and the metadata variable 
# "genotype" contains "Control" and "KO". Thus, 4 colors are needed for the
# following cell groups: Female_Control, Female_KO, Male_Control, Male_KO
# plot_colors <- c("#ef4444", "#06b6d4", "#f59e0b", "#0d9488")

#------------------------------------------------------------------------------#
#----                              Main script                             ----#
#------------------------------------------------------------------------------#

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

# Plot cells in UMAP space, colored by metadata variable(s)
# Note: to plot cells in t-SNE space, replace instances of 'umap'/'UMAP' with 
# 'tsne'/'tSNE'. To plot cells in PCA space, replace instances of 'umap'/'UMAP'
# with 'pca'/'PC'.

# Add new column to Seurat object meta.data for cell grouping
grouping_col <- paste(cell_grouping, collapse = "_")
if (length(cell_grouping) > 1) {
    so@meta.data[[grouping_col]] <- apply(
        so@meta.data[, cell_grouping], 1, paste, collapse = "_"
    )
}

# Get plot data (target expression and coordinates)
plot_data <- bind_cols(
    as.data.frame(so@reductions$umap@cell.embeddings),
    as.data.frame(so@meta.data[[grouping_col]])
)
colnames(plot_data)[ncol(plot_data)] <- grouping_col

# Create cell scatter plot
cell_scatter_plot <- plot_ly(plot_data) %>%
    add_trace(
        x = plot_data$UMAP_1,
        y = plot_data$UMAP_2,
        type = "scattergl",
        mode = "markers",
        color = plot_data[[grouping_col]],
        colors = plot_colors,
        marker = list(size = 3),
        legendgroup = plot_data[[grouping_col]],
        showlegend = TRUE
    ) %>%
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
    title = paste0(
        "Cell scatter plot, colored by ", 
        paste(cell_grouping, collapse = ", ")
    ),
    selfcontained = TRUE
)

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)