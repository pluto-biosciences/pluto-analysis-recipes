################################################################################
#####                   Differential expression analysis                   #####
################################################################################

# Script to compare genome-wide expression profiles of two groups of cells and
# create a volcano plot visualization using data in a Pluto Bio experiment.

# Note that this script is designed for scRNA-seq data.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

# For more advanced options with running Seurat `FindMarkers()`, see the Seurat
# usage docs: https://satijalab.org/seurat/reference/findmarkers

# For more advanced plot customization, scroll down to the Main script and edit
# parameters for the  plot_ly function. For more information, see the Plotly
# usage docs: https://plotly.com/r/line-and-scatter/

################################################################################
#####                        User-defined parameters                       #####
################################################################################

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX104142"

# List the Seurat object meta.data columns you want to use to group cells by.
# You can view these by running `head(so@meta.data)`.
metadata_cols <- c("cellkb_annotations", "disease_state")

# Define the target group of cells you want to use for the analysis.
# In target_group, separate groups with an underscore.
# target_group_name is the display name for the target group.
target_group <- "B cell 1_Chronic Colitis"
target_group_name <- "B cell 1 - Chronic Colitis"

# Define the reference group of cells you want to use for the analysis.
# In reference_group, separate groups with an underscore.
# reference_group_name is the display name for the reference group.
reference_group <- "B cell 1_Healthy"
reference_group_name <- "B cell 1 - Healthy"

# Statistics options
# Test to use for Seurat `FindMarkers()` function:
stat_test <- "wilcox"
# Only test genes that are detected in the minimum fraction of cells:
min_pct <- 0.01

# Define a name for the analysis in Pluto.
analysis_name <- "Differential expression analysis"

# Define a file path for the analysis plot (volcano plot).
display_file_path <- paste0(
    experiment_id, "_deg_plot_",
    target_group_name, "_vs_", reference_group_name,
    ".html"
)

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(
    experiment_id, "_deg_results_",
    target_group_name, "_vs_", reference_group_name,
    ".csv"
)

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Differential expression analysis was performed with the Seurat R package,",
    " comparing the groups: ", target_group_name, " vs ", reference_group_name,
    ". Differentially expressed genes between ", target_group_name, " vs ",
    reference_group_name, " were identified using the ", stat_test, " test. ",
    "Genes were only tested if they were dtected in a minimum fraction of ",
    min_pct, " cells in either of the two groups."
)

# Define plot colors.
sig_increased_color <- "#EF4444"
sig_decreased_color <- "#06B6D4"
not_sig_color <- "gray"

# Define plot significance thresholds.
pval_threshold <- 0.01
lower_fc_threshold <- 0.8 # Fold change; will be converted to log2
upper_fc_threshold <- 1.2 # Fold change; will be converted to log2

################################################################################
#####                          Helper function(s)                          #####
################################################################################

get_plot_data <- function(
    deg_table, pval_threshold, lower_fc_threshold, upper_fc_threshold) {
    lower_log2fc_threshold <- log2(lower_fc_threshold)
    upper_log2fc_threshold <- log2(upper_fc_threshold)
    plot_data <- deg_table %>%
        mutate(
            # Add a new column for the negative log10 adjusted p-value
            neg_log10_pvalue = -log10(Adj_P_Value),
            # Replace Inf values with 310 for plotting purposes
            neg_log10_pvalue = if_else(
                is.infinite(neg_log10_pvalue), 310, neg_log10_pvalue
            ),
            # Add a new column for the significance category of each point
            sig_category = case_when(
                Adj_P_Value <= pval_threshold &
                    Log2_Fold_Change < lower_log2fc_threshold ~ "sig_decreased",
                Adj_P_Value <= pval_threshold &
                    Log2_Fold_Change > upper_log2fc_threshold ~ "sig_increased",
                TRUE ~ "not_sig"
            )
        )
    return(plot_data)
}

################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(Seurat)
library(tidyverse)
library(plotly)
library(htmlwidgets)

#------------------------------------------------------------------------------#
#----                             Prepare data                             ----#
#------------------------------------------------------------------------------#
# Log into Pluto
pluto_login(api_token)

# Read in experiment Seurat object
final_obj <- pluto_read_seurat_object(experiment_id, "final")
so <- final_obj$obj

#------------------------------------------------------------------------------#
#----                        Differential expression                       ----#
#------------------------------------------------------------------------------#
# Add new grouping column to Seurat object meta.data
message("Adding in new deg_comp column to Seurat meta.data...")
message("Metadata cols: ", paste(metadata_cols, collapse = ", "))
if (length(metadata_cols) > 1) {
    so@meta.data$deg_comp <- apply(
        so@meta.data[, metadata_cols],
        1,
        paste,
        collapse = "_"
    )
} else {
    so@meta.data$deg_comp <- so@meta.data[[metadata_cols[[1]]]]
}

# Set active identity to new deg_comp column
message("Setting active Seurat identity to deg_comp...")
Idents(object = so) <- so@meta.data$deg_comp


# Find differentially expressed genes
message("Finding DEGs for ", target_group, " vs ", reference_group, "...")
message("Statistical test: ", stat_test)
deg_table <- FindMarkers(
    so,
    assay = "RNA",
    ident.1 = target_group,
    ident.2 = reference_group,
    log2fc.threshold = 0,
    test.use = stat_test,
    min.pct = min_pct,
    only.pos = FALSE
)

# Format DEG table
message("Formatting DEG table...")
deg_table <- tibble::rownames_to_column(deg_table, "Gene_Symbol")
col_order <- c(
    "Gene_Symbol",
    "avg_log2FC",
    "p_val",
    "p_val_adj",
    "pct.1",
    "pct.2"
)
deg_table <- deg_table[, col_order] %>%
    dplyr::rename(
        "Log2_Fold_Change" = "avg_log2FC",
        "P_Value" = "p_val",
        "Adj_P_Value" = "p_val_adj",
        "pct_1" = "pct.1",
        "pct_2" = "pct.2"
    )

# Save DEG results table as .csv file
write.csv(deg_table, results_file_path, row.names = FALSE, na = "")

#------------------------------------------------------------------------------#
#----                         Create a volcano plot                        ----#
#------------------------------------------------------------------------------#
# Create a volcano plot, Pluto style!
plot_data <- get_plot_data(
    deg_table, pval_threshold, lower_fc_threshold, upper_fc_threshold
)

# Define custom color palette
color_palette <- c(
    "sig_increased" = sig_increased_color,
    "sig_decreased" = sig_decreased_color,
    "not_sig" = not_sig_color
)

# Create the volcano plot
volcano_plot <- plot_ly(
    data = plot_data,
    x = ~Log2_Fold_Change,
    y = ~neg_log10_pvalue,
    type = "scatter",
    mode = "markers",
    text = ~Gene_Symbol, # Hover text with Gene Symbol
    hoverinfo = "text", # Show only gene symbol on hover
    color = ~ factor(sig_category), # Color by significance category
    colors = color_palette, # Custom color palette
    marker = list(
        size = 9, # Point size
        opacity = 0.75 # Point opacity
    )
) %>%
    layout(
        xaxis = list(title = "log2 fold change", showgrid = FALSE),
        yaxis = list(title = "-log10(adjusted p-value)", showgrid = FALSE),
        showlegend = FALSE
    )

# Save volcano plot
saveWidget(
    volcano_plot,
    file = display_file_path,
    title = paste0(target_group_name, " vs ", reference_group_name),
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
