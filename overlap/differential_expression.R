################################################################################
#####           Overlap analysis - differentially expressed genes          #####
################################################################################

# Script to create a Venn diagram of overlapping genes from multiple
# differential expression analyses in Pluto Bio.

# This example script creates a Venn diagram for overlapping upregulated genes.
# To create a Venn diagram for overlapping downregulated genes, adjust the
# filtering thresholds and direction.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

# For more advanced plot customization, scroll down to the Main script and edit
# parameters for the venn.diagram function. For more information, see the usage
# docs: https://cran.r-project.org/web/packages/VennDiagram/VennDiagram.pdf

################################################################################
#####                        User-defined parameters                       #####
################################################################################

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX276186"

# Define a file path for the analysis plot (Venn diagram).
display_file_path <- paste0(experiment_id, "_venn_diagram.png")

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(experiment_id, "_overlap_results.csv")

# Define a name for the analysis in Pluto.
analysis_name <- "Overlap analysis: Upregulated genes"

# Include methods to describe your analysis.
plot_methods <- paste(
    "Overlapping upregulated genes,",
    "plotted using the `R` `VennDiagram` package."
)

# List the differential expression analyses you want to include.
# You will need the plot ID for each analysis.
# The plot ID for an analysis is listed under the Methods section for that
# analysis in Pluto.
plot_ids <- c(
    "3b37562f-e0e2-4dbe-8aed-e3ebda2076df",
    "f26be9c4-44b2-464e-a33c-66756f28b482",
    "c9324c25-1a7e-44fa-8c41-29db892f2af9"
)

# List category names that correspond to each of the above analyses.
# Note: order matters! The order of category_names should match the order of
# plot_ids.
category_names <- c(
    "80epi",
    "12ss",
    "24hpf"
)

# Define plot colors (one color per differential expression analysis).
# Note: the length of plot_colors must match the length of plot_ids and
# category_names.
plot_colors <- c(
    "#06141B",
    "#AD05B1",
    "#4B1590"
)

# Define filtering parameters to get significant genes.
adj_pvalue_threshold <- 0.05
fold_change_threshold <- 1.5
filter_direction <- "up" # use "down" for downregulated genes

################################################################################
#####                          Helper function(s)                          #####
################################################################################

filter_sig_genes <- function(
    df, adj_pvalue_threshold, fold_change_threshold, filter_direction) {
    if (filter_direction == "up") {
        df %>%
            filter(
                Adj_P_Value <= adj_pvalue_threshold,
                Log2_Fold_Change > log2(fold_change_threshold)
            )
    } else if (filter_direction == "down") {
        df %>%
            filter(
                Adj_P_Value <= adj_pvalue_threshold,
                Log2_Fold_Change < log2(fold_change_threshold)
            )
    } else {
        stop("Invalid parameter! filter_direction must be 'up' or 'down'.")
    }
}

strip_trailing_digits <- function(x, category_names) {
    pattern <- paste0(
        "(?<=",
        paste(category_names, collapse = "|"),
        ")[0-9]+$"
    )
    gsub(pattern, "", x, perl = TRUE)
}


################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(tidyverse)
library(VennDiagram)
library(gplots)

# Log into Pluto
pluto_login(api_token)

# Read in analysis results
results_list <- lapply(plot_ids, function(plot_id) {
    tryCatch(
        {
            pluto_read_results(experiment_id = experiment_id, plot_id = plot_id)
        },
        error = function(e) {
            stop(
                "Failed to read results for experiment ID: ", experiment_id,
                " and plot ID: ", plot_id
            )
        }
    )
})

# Filter results for significant genes
filtered_results <- lapply(results_list, function(res) {
    filter_sig_genes(
        res, adj_pvalue_threshold, fold_change_threshold, filter_direction
    )
})

# Get lists of significant genes
gene_lists <- lapply(filtered_results, function(res) {
    res$Gene_Symbol
})

# Plot Venn diagram for significant genes
venn_plot <- venn.diagram(
    x = gene_lists,
    category.names = category_names,
    fill = plot_colors,
    lwd = 2,
    lty = "blank",
    cex = 1,
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    cat.default.pos = "outer",
    filename = NULL
)
png(
    display_file_path,
    width = 1200,
    height = 1200,
    res = 300
)
grid.draw(venn_plot)
dev.off()

# Get Venn diagram results in tabular format
names(gene_lists) <- category_names
gene_sets <- venn(gene_lists, show.plot = FALSE)
overlap_results <- unlist(attributes(gene_sets)$intersections)
overlap_results <- data.frame(
    gene_symbol = overlap_results,
    set = names(overlap_results)
)
overlap_results$set <- strip_trailing_digits(
    overlap_results$set,
    category_names
)
write.csv(overlap_results, results_file_path, row.names = FALSE)

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    results_file_path = results_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)
