################################################################################
#####              Marker expression analysis - single marker              #####
################################################################################

# Script to create visualizations for expression of a single marker gene for
# single cell RNA-seq (scRNA-seq) data in Pluto Bio.

# This example script illustrates how to create a cell scatter plot in UMAP
# space and a violin plot for a given marker gene. The cell scatter plot is
# equivalent to the Feature plot in Seurat (https://satijalab.org/seurat/), and
# can utilize t-SNE or PCA coordinates in place of UMAP coordinates.

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


################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(tidyverse)
library(VennDiagram)

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

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)


############################## WHEEEEEEEEEEEEEEEEEEEEEE###################################################
library(pluto)
library(tidyverse)
library(plotly)
library(htmlwidgets)
library(ggplot2)
library(Seurat)

# Log into Pluto
# pluto_login('<YOUR_API_TOKEN>')
pluto_login("489738dfae03331601d9e44a593a2c252d760126")

# Get data from Pluto
experiment_id <- "PLX220431"

# exp_obj <- pluto_read_seurat_object(
#  experiment_id = experiment_id,
#  seurat_type = 'final'
# )
so <- readRDS("seurat_object.rds")

# Downsample for plotting
# so.downsample <- so[, sample(colnames(so), size = 100000, replace = FALSE)]

# Define genes of interest for plotting
targets <- read.csv("postera_target_genes_updated.csv")
targets <- targets$gene_symbol

# Plot target expression in UMAP space
for (target in targets) {
    if (target %in% rownames(so@assays$RNA@data)) {
        # Create UMAP plot, colored by gene expression
        umap_data <- bind_cols(
            as.data.frame(so@reductions$umap@cell.embeddings),
            as.data.frame(so@assays$RNA@data[`target`, ])
        )
        colnames(umap_data)[ncol(umap_data)] <- target

        umap <- plot_ly(umap_data) %>%
            add_trace(
                x = umap_data$UMAP_1,
                y = umap_data$UMAP_2,
                type = "scattergl",
                mode = "markers",
                color = umap_data[[target]],
                colors = c("#cdcfd1", "#AD05B1"),
                marker = list(size = 3),
                legendgroup = umap_data[[target]],
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

        # Save UMAP plot
        saveWidget(
            umap,
            file = paste0("plots/", experiment_id, "_umap_", target, ".html"),
            title = paste0("UMAP, colored by ", target, " expression"),
            selfcontained = TRUE
        )

        # Push results back to Pluto
        pluto_add_experiment_plot(
            experiment_id = experiment_id,
            display_file_path = paste0("plots/", experiment_id, "_umap_", target, ".html"),
            analysis_name = paste0("UMAP, colored by ", target, " expression"),
            plot_methods = paste0("Cells in UMAP space, colored by ", target, " expression.")
        )
    } else {
        print(paste0(target, " is not present in the Seurat object assay data."))
    }
}

# Plot UMAP colored by default clustering resolution 0.2
umap_data <- bind_cols(
    as.data.frame(so@reductions$umap@cell.embeddings),
    as.data.frame(so@meta.data)
)

custom_colors <- c(
    "#6366f1", "#a5b4fc", "#9333ea", "#c084fc", "#38bdf8", "#0369a1",
    "#8b5cf6", "#c4b5fd", "#4f46e5", "#818cf8", "#a855f7", "#d8b4fe",
    "#0284c7", "#7dd3fc", "#7c3aed", "#a78bfa", "#22d3ee", "#14b8a6",
    "#fcd34d", "#f97316", "#ef4444", "#06b6d4", "#fbbf24", "#0d9488",
    "#ea580c", "#67e8f9", "#2dd4bf", "#fde68a", "#fb923c", "#f87171",
    "#5eead4", "#f59e0b", "#0891b2", "#c2410c", "#312e81", "#fca5a5",
    "#0c4a6e", "#0e7490", "#0f766e", "#059669", "#65a30d", "#84cc16",
    "#f59e0b", "#fdba74", "#fecdd3", "#f5d0fe", "#7f1d1d", "#0284c7",
    "#b45309", "#ca8a04", "#fde047", "#fef9c3", "#ccfbf1", "#7dd3fc",
    "#0ea5e9", "#0369a1", "#1e40af"
)

umap <- plot_ly(umap_data) %>%
    add_trace(
        x = umap_data$UMAP_1,
        y = umap_data$UMAP_2,
        type = "scattergl",
        mode = "markers",
        color = umap_data$res_0.1,
        colors = custom_colors,
        marker = list(size = 3),
        legendgroup = umap_data$res_0.1,
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
    style(hoverinfo = "none")

# Save UMAP plot
saveWidget(
    umap,
    file = paste0("plots/", experiment_id, "_umap_res_0_1.html"),
    title = paste0("UMAP, colored by clustering resolution 0.1"),
    selfcontained = TRUE
)

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = paste0("plots/", experiment_id, "_umap_res_0_1.html"),
    analysis_name = paste0("UMAP, colored by clustering resolution 0.1"),
    plot_methods = paste0("Cells in UMAP space, colored by clustering resolution 0.1.")
)

# Plot UMAP colored by major class
major_class_colors <- c("#8b5cf6", "#0284c7", "#14b8a6", "#fcd34d", "#ea580c")
umap <- plot_ly(umap_data) %>%
    add_trace(
        x = umap_data$UMAP_1,
        y = umap_data$UMAP_2,
        type = "scattergl",
        mode = "markers",
        color = umap_data$major_class,
        colors = major_class_colors,
        marker = list(
            size = 3
        ),
        legendgroup = umap_data$major_class,
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
    style(hoverinfo = "none")

# Save UMAP plot
saveWidget(
    umap,
    file = paste0("plots/", experiment_id, "_umap_major_class.html"),
    title = paste0("UMAP, colored by major class"),
    selfcontained = TRUE
)

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = paste0("plots/", experiment_id, "_umap_major_class.html"),
    analysis_name = paste0("UMAP, colored by major class"),
    plot_methods = paste0("Cells in UMAP space, colored by major class.")
)

# Plot UMAP colored by anatomical location
umap <- plot_ly(umap_data) %>%
    add_trace(
        x = umap_data$UMAP_1,
        y = umap_data$UMAP_2,
        type = "scattergl",
        mode = "markers",
        color = umap_data$anatomical_location,
        colors = custom_colors,
        marker = list(
            size = 3
        ),
        legendgroup = umap_data$anatomical_location,
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
    style(hoverinfo = "none")

# Save UMAP plot
saveWidget(
    umap,
    file = paste0("plots/", experiment_id, "_umap_anatomical_location.html"),
    title = paste0("UMAP, colored by anatomical location"),
    selfcontained = TRUE
)

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = paste0("plots/", experiment_id, "_umap_anatomical_location.html"),
    analysis_name = paste0("UMAP, colored by anatomical location"),
    plot_methods = paste0("Cells in UMAP space, colored by anatomical location.")
)

# Plot target expression as violin plots
for (target in targets) {
    if (target %in% rownames(so@assays$RNA@data)) {
        # Create violin plot, grouped by major class
        vln_plt <- so@meta.data %>%
            mutate(
                expression = so@assays$RNA@data[`target`, ],
                expression = expression + rnorm(nrow(.)) / 200
            ) %>%
            ggplot(aes(x = major_class, y = expression, fill = major_class)) +
            geom_violin(draw_quantiles = c(0.5), scale = "width", trim = TRUE) +
            theme_bw() +
            scale_fill_manual(values = major_class_colors) +
            scale_x_discrete() +
            scale_y_continuous(
                name = "Log-normalized expression",
                labels = scales::comma,
                limits = c(0, 0.25)
            ) +
            labs(title = paste0(target, " expression by major class")) +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                legend.position = "none"
            )

        # Save violin plot
        ggsave(
            paste0("plots/", experiment_id, "_violin_plot_by_major_class_", target, "_zoom.png"),
            height = 6, width = 6
        )

        # Push results back to Pluto
        pluto_add_experiment_plot(
            experiment_id = experiment_id,
            display_file_path = paste0("plots/", experiment_id, "_violin_plot_by_major_class_", target, "_zoom.png"),
            analysis_name = paste0("Violin plot, ", target, " expression grouped by major class (zoom)"),
            plot_methods = paste0("Violin plot, ", target, " expression grouped by major class. Violins are scaled by width and plotted with additional noise to avoid highlighting distinct values at the bottom end of the scale (done by default with Seurat). Additionally, violins are trimmed to avoid the violin extending into the negative scale (as this does not make sense for expression values).")
        )
    } else {
        print(paste0(target, " is not present in the Seurat object assay data."))
    }
}

# Plot all targets as dot plot
dot_plt <- DotPlot(
    so,
    features = unique(targets[targets != "IL8"]),
    group.by = "major_class"
) +
    labs(x = "Target genes", y = "Major class") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    )

ggsave(
    paste0("plots/", experiment_id, "_dot_plot_all_targets.png"),
    height = 6, width = 12
)

# Push dot plot back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = paste0("plots/", experiment_id, "_dot_plot_all_targets.png"),
    analysis_name = paste0("Dot plot, all targets grouped by major class"),
    plot_methods = paste0("Dot plot, highlighting target gene expression across major classes. The size of the dot encodes the percentage of cells within a class, while the color encodes the average expression level across all cells within a class.")
)

# Make table for target genes
counts <- so@assays$RNA@counts
major_class <- unique(so@meta.data$major_class)
counts_list <- list()

for (class in major_class) {
    class_cells <- which(so@meta.data$major_class == class)
    class_counts <- counts[, class_cells]
    ncells <- ncol(class_counts)

    counts_df <- data.frame(
        target_gene = character(),
        total_cells = numeric(),
        zero_expr_count = numeric(),
        nonzero_expr_count = numeric(),
        ratio_expr = numeric(),
        pct_expr = numeric()
    )

    for (target in targets) {
        if (target %in% row.names(class_counts)) {
            nonzero_count <- sum(class_counts[target, ] > 0)
            zero_count <- ncells - nonzero_count
            ratio <- nonzero_count / ncells
            pct <- (nonzero_count / ncells) * 100
            counts_df <- rbind(
                counts_df,
                data.frame(
                    target_gene = target,
                    total_cells = ncells,
                    zero_expr_count = zero_count,
                    nonzero_expr_count = nonzero_count,
                    ratio_expr = ratio,
                    pct_expr = pct
                )
            )
        } else {
            print(paste0(target, " is not present in the Seurat object assay data."))
        }
    }

    counts_list[[class]] <- counts_df
}

combined_counts_df <- bind_rows(counts_list, .id = "major_class")
combined_counts_df$pct_expr <- paste0(format(round(combined_counts_df$pct_expr, 2), nsmall = 2), "%")

write.csv(combined_counts_df, "gene_expression_cell_counts_by_major_class.csv", row.names = FALSE)
