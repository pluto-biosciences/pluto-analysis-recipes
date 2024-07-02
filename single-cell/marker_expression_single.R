################################################################################
#####              Marker expression analysis - single marker              #####
################################################################################

# Script to create cell scatter plot and violin plot visualizations for
# expression of a single marker gene for single cell RNA-seq (scRNA-seq) data in
# Pluto Bio.

################################################################################
#####                           Cell scatter plot                          #####
################################################################################

# This example script illustrates how to create a cell scatter plot in UMAP
# space for a given marker gene. The cell scatter plot is equivalent to the
# Feature plot in Seurat (https://satijalab.org/seurat/), and can utilize t-SNE
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

################################################################################
#####                              Violin plot                             #####
################################################################################

# This example script illustrates how to create a violin plot, grouping cells
# by the timepoint variable, which is available in this experiment's sample
# data/Seurat meta.data. You can group your data by any variable available in
# your Seurat meta.data. To see available variables, run:
# head(so@meta.data)
# once you read in your Seurat object in the Main script.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

# For more advanced plot customization, scroll down to the Main script (Violin
# plot) and edit parameters for the  ggplot function. For more information, see
# the usage docs: https://ggplot2.tidyverse.org/index.html

#------------------------------------------------------------------------------#
#----                        User-defined parameters                       ----#
#------------------------------------------------------------------------------#

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX141180"

# Define your target marker gene.
target <- "olig2"

# Define a file path for the analysis plot (violin plot).
display_file_path <- paste0(experiment_id, "_", target, "_violin.png")

# Define a name for the analysis in Pluto.
analysis_name <- paste0("Marker expression analysis: ", target)

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Cell scatter plot using UMAP coordinates, colored by ",
    target, " expression."
)

# Define your grouping variable.
# Remember, you can access available grouping variables by running:
# head(so@meta.data)
cell_grouping <- "timepoint"

# Colors for violin plot groups.
# The length of plot_colors must match the number of categories contained
# in cell_grouping.
plot_colors <- c("#AD05B1", "#f59e0b", "#0d9488")

# If you are grouping cells by one of your clustering resolutions, you can
# access the same color palette you are using in the Pluto app by loading in
# your custom color palettes and selecting the palette that corresponds to your
# resolution of interest after loading in exp_obj (see Main script).
# custom_palettes <- exp_obj$colors
# plot_colors <- custom_palettes$default_clusters_res_0_1
# For more info, see:
# https://pluto-biosciences.github.io/pluto-sdk-R/articles/scrnaseq_recipes.html

#------------------------------------------------------------------------------#
#----                              Main script                             ----#
#------------------------------------------------------------------------------#

# Load required libraries
library(pluto)
library(tidyverse)
library(ggplot2)

# Log into Pluto
pluto_login(api_token)

# Get Seurat object
exp_obj <- pluto_read_seurat_object(
    experiment_id = experiment_id,
    seurat_type = "final"
)
so <- exp_obj$obj

# Create violin plot
vln_plot <- so@meta.data %>%
    mutate(
        expression = so@assays$RNA@data[`target`, ],
        expression = expression + rnorm(nrow(.)) / 200
    ) %>%
    ggplot(aes(
        x = !!sym(cell_grouping),
        y = expression,
        fill = !!sym(cell_grouping)
    )) +
    geom_violin(draw_quantiles = c(0.5), scale = "width", trim = TRUE) +
    theme_bw() +
    scale_fill_manual(values = plot_colors) +
    scale_x_discrete() +
    scale_y_continuous(
        name = "Log-normalized expression",
        labels = scales::comma
    ) +
    labs(title = paste0(target, " expression by ", cell_grouping)) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none"
    )

# Save violin plot
ggsave(
    display_file_path,
    height = 6, width = 6
)

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)
