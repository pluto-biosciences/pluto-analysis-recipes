# Comprehensive Script to Create Venn Diagrams for Gene Expression Analysis
# in Pluto Bio

# Load required libraries
library(pluto)
library(tidyverse)
library(VennDiagram)

# Helper functions to ensure secure and robust operations
safe_pluto_login <- function() {
    api_token <- Sys.getenv("PLUTO_API_TOKEN")
    if (nchar(api_token) == 0) stop("API token is missing")
    pluto_login(api_token)
}

safe_pluto_read_results <- function(experiment_id, plot_id) {
    tryCatch(
        {
            pluto_read_results(experiment_id = experiment_id, plot_id = plot_id)
        },
        error = function(e) {
            message("Failed to read results for experiment ID: ", experiment_id, " and plot ID: ", plot_id)
            NULL # return NULL on error
        }
    )
}

filter_significant_genes <- function(df, adj_pvalue_threshold, log2fc_threshold) {
    df %>%
        filter(Adj_P_Value < adj_pvalue_threshold, Log2_Fold_Change >= log2fc_threshold)
}

# Log into Pluto safely
safe_pluto_login()

# Define experiment IDs and plot IDs
experiment_id <- "PLX276186"
plot_ids <- c("3b37562f-e0e2-4dbe-8aed-e3ebda2076df", "f26be9c4-44b2-464e-a33c-66756f28b482", "c9324c25-1a7e-44fa-8c41-29db892f2af9")

# Read results safely and filter for significant genes
adj_pvalue_threshold <- 0.05
log2fc_threshold <- log2(1.5)

results_list <- lapply(plot_ids, function(plot_id) {
    safe_pluto_read_results(experiment_id, plot_id)
})

# Extract upregulated genes
upreg_results <- lapply(results_list, function(res) {
    filter_significant_genes(res, adj_pvalue_threshold, log2fc_threshold)
})

# Generate lists of upregulated genes
upreg_gene_lists <- lapply(upreg_results, function(res) {
    res$Gene_Symbol
})

# Plot Venn diagram for upregulated genes
venn_upreg <- venn.diagram(
    x = upreg_gene_lists,
    category.names = c("80epi", "12ss", "24hpf"),
    fill = c("#06141B", "#AD05B1", "#4B1590"),
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
    paste0("plots/", experiment_id, "_venn_diagram_upreg.png"),
    width = 1200,
    height = 1200,
    res = 300
)
grid.draw(venn_upreg)
dev.off()

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = paste0("plots/", experiment_id, "_venn_diagram_upreg.png"),
    analysis_name = "Overlap analysis: Upregulated genes",
    plot_methods = "Overlapping upregulated genes, plotted using the `R` `VennDiagram` package."
)

# Note: Downregulated genes analysis would follow a similar pattern with appropriate filters and plotting.
