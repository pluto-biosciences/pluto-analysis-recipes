################################################################################
#####                     Over-representation analysis                     #####
################################################################################

# Script to run over-representation analysis (ORA) and create a dot plot of the
# top 20 pathways or gene sets that are most significantly over-represented in a
# given list of differentially expressed genes from Pluto Bio.

# This example script runs ORA on a list of genes identified as significantly
# upregulated in differential expression analysis. To run ORA on downregulated
# genes or on all differential genes, adjust the filtering thresholds and
# direction defined under User-defined parameters.

# Additionally, this script highlights two approaches to running ORA. One
# approach shows how to use a gene set collection availabe in MSigDB. The other
# approach shows how to use a custom gene set collection.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

# For more advanced plot customization, scroll down to the Main script and edit
# parameters for the ggplot function. For more information, see the usage docs:
# https://ggplot2.tidyverse.org/

################################################################################
#####                      MSigDB gene set collection                      #####
################################################################################

#------------------------------------------------------------------------------#
#----                        User-defined parameters                       ----#
#------------------------------------------------------------------------------#

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX276186"

# Define a file path for the analysis plot (dot plot).
display_file_path <- paste0(experiment_id, "_ORA_dot_plot.png")

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(experiment_id, "_ORA_results.csv")

# Define a name for the analysis in Pluto.
analysis_name <- "Over-representation analysis: Upregulated genes"

# Speficy your experiment's organism.
# Acceptable organisms include: human, mouse, rat, zebrafish, fruit fly, xlaevis
organism <- "zebrafish"

# Specify the collection shortname you would like to use from MSigDB.
# Note: to access collection shortnames programmatically, run
# unique(all_gene_sets$Collection_Shortname) after the prepare gene sets step.
collection_shortname <- "c5__go_bp"

# Specify the differential expression analysis you want to use for your gene
# list. This list will be filtered using adjusted p-value and fold change
# thresholds.
plot_id <- "f26be9c4-44b2-464e-a33c-66756f28b482"

# Define the adjusted p-value threshold to filter your gene list by.
adj_pvalue_threshold <- 0.01

# Define the fold change direction to filter your gene list by.
fold_change_direction <- "up" # up, down, or both

# Define the fold change threshold to filter your gene list by.
fold_change_threshold <- 1.2

# Define the minimum number of genes that a gene set must contain to be
# included in the analysis.
min_gene_set_size <- 5

# Define the maximum number of genes that a gene set must contain to be
# included in the analysis.
max_gene_set_size <- 1000

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Over-representation analysis was performed using the `clusterProfiler` R",
    " package. Passing genes were defined as those with an adjusted p-value of",
    " less than or equal to ", adj_pvalue_threshold, ", and a log fold-change",
    " greater than ", fold_change_threshold, "."
)

# Define plot colors for dots (colored by adjusted p-value).
low_color <- "blue"
high_color <- "red"

#------------------------------------------------------------------------------#
#----                          Helper function(s)                          ----#
#------------------------------------------------------------------------------#

get_all_gene_sets <- function(organism) {
    all_gene_sets <- msigdbr(species = organism)
    all_gene_sets$gs_subcat <- gsub("\\:", "_", all_gene_sets$gs_subcat)
    all_gene_sets$Collection_Shortname <- tolower(paste0(
        all_gene_sets$gs_cat,
        ifelse(all_gene_sets$gs_subcat == "", "", "__"),
        all_gene_sets$gs_subcat
    ))
    names(all_gene_sets)[3] <- "Gene_Set_Shortname"
    return(all_gene_sets)
}

check_target_gene_set <- function(all_gene_sets, collection_shortname) {
    unique_collections <- unique(all_gene_sets$Collection_Shortname)
    if (!collection_shortname %in% unique_collections) {
        stop("Unsupported collection_shortname: ", collection_shortname,
            call. = FALSE
        )
    }
}

get_target_gene_sets <- function(all_gene_sets, collection_shortname) {
    target_gene_sets <- unique(
        all_gene_sets[
            all_gene_sets$Collection_Shortname == collection_shortname,
        ]
    )
    target_gene_sets <- target_gene_sets[
        , c("Gene_Set_Shortname", "gene_symbol")
    ]
    return(target_gene_sets)
}

get_gene_list <- function(
    df, adj_pvalue_threshold, fold_change_threshold, fold_change_direction) {
    if (fold_change_direction == "up") {
        filtered_df <- df %>%
            filter(
                Adj_P_Value <= adj_pvalue_threshold,
                Log2_Fold_Change > log2(fold_change_threshold)
            )
    } else if (fold_change_direction == "down") {
        filtered_df <- df %>%
            filter(
                Adj_P_Value <= adj_pvalue_threshold,
                Log2_Fold_Change < log2(fold_change_threshold)
            )
    } else if (fold_change_direction == "both") {
        filtered_df <- df %>%
            filter(
                Adj_P_Value <= adj_pvalue_threshold,
                Log2_Fold_Change < log2(fold_change_threshold) |
                    Log2_Fold_Change < -log2(fold_change_threshold)
            )
    }
    gene_list <- filtered_df$Gene_Symbol
    return(unique(gene_list))
}

format_ora_table <- function(ora_obj) {
    ora_table <- data.frame(ora_obj)
    ora_table <- ora_table[, -grep("Description", colnames(ora_table))]
    ora_table <- ora_table[order(ora_table$p.adjust, decreasing = FALSE), ]
    return(ora_table)
}

get_plot_data <- function(ora_table) {
    plot_data <- ora_table %>% slice_head(n = 20)
    plot_data <- plot_data %>%
        mutate(
            GeneRatio = sapply(GeneRatio, function(x) {
                parts <- strsplit(x, "/")[[1]]
                if (length(parts) == 2) {
                    as.numeric(parts[1]) / as.numeric(parts[2])
                } else {
                    NA # Handle cases where the format might be incorrect
                }
            })
        )
    # Reorder ID by GeneRatio
    plot_data$ID <- factor(
        plot_data$ID,
        levels = plot_data$ID[order(plot_data$GeneRatio)]
    )
    return(plot_data)
}

#------------------------------------------------------------------------------#
#----                              Main script                             ----#
#------------------------------------------------------------------------------#

# Load required libraries
library(pluto)
library(clusterProfiler)
library(msigdbr)
library(tidyverse)

# Log into Pluto
pluto_login(api_token)

# Read in analysis results
deg_table <- pluto_read_results(
    experiment_id = experiment_id,
    plot_id = plot_id
)

# Prepare gene sets
all_gene_sets <- get_all_gene_sets(organism)
# unique(all_gene_sets$Collection_Shortname)

# Check that provided collection_shortname is supported
check_target_gene_set(all_gene_sets, collection_shortname)

# Get target gene sets
target_gene_sets <- get_target_gene_sets(all_gene_sets, collection_shortname)

# Generate gene list from differential expression table
gene_list <- get_gene_list(
    deg_table, adj_pvalue_threshold,
    fold_change_threshold, fold_change_direction
)

# Run ORA
ora_obj <- enricher(
    gene_list,
    TERM2GENE = target_gene_sets,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    qvalueCutoff = 1,
    pvalueCutoff = 1,
    universe = deg_table$Gene_Symbol,
    pAdjustMethod = "BH"
)

# Format ORA table
ora_table <- format_ora_table(ora_obj)

# Save ORA results
write.csv(ora_table, results_file_path, row.names = FALSE, na = "")

# Generate dot plot for top 20 pathways or gene sets
plot_data <- get_plot_data(ora_table)

plt <- ggplot(
    plot_data,
    aes(x = GeneRatio, y = ID, color = p.adjust, size = Count)
) +
    geom_point() +
    scale_color_gradient(low = low_color, high = high_color) +
    scale_size_continuous(range = c(2, 10)) +
    labs(
        x = "Gene Ratio",
        y = "Gene Set",
        color = "Adj. p-value",
        size = "Count"
    ) +
    theme_minimal() +
    theme(
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "in"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)
    )

ggsave(
    display_file_path,
    plot = plt,
    width = 12, height = 9, units = "in", dpi = 300
)

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)

################################################################################
#####                      Custom gene set collection                      #####
################################################################################

#------------------------------------------------------------------------------#
#----                        User-defined parameters                       ----#
#------------------------------------------------------------------------------#

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX276186"

# Define a file path for the analysis plot (dot plot).
display_file_path <- paste0(experiment_id, "_ORA_dot_plot.png")

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(experiment_id, "_ORA_results.csv")

# Define a name for the analysis in Pluto.
analysis_name <- "Over-representation analysis: Upregulated genes"

# Specify the path to your custom gene set collection CSV file.
target_gene_set_file <- "custom_gene_set_collection.csv"
# Note: target_gene_set_file must be formatted with a "Gene_Set_Shortname"
# column and a "gene_symbol" columm, such as in this example:
# Gene_Set_Shortname,gene_symbol
# GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS,AASDHPPT
# GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS,ALDH1L1
# GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS,ALDH1L2
# GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS,MTHFD1
# GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS,MTHFD1L
# GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS,MTHFD2L
# GOBP_2FE_2S_CLUSTER_ASSEMBLY,BOLA2
# GOBP_2FE_2S_CLUSTER_ASSEMBLY,BOLA2B
# GOBP_2FE_2S_CLUSTER_ASSEMBLY,GLRX3
# GOBP_2FE_2S_CLUSTER_ASSEMBLY,GLRX5
# GOBP_2FE_2S_CLUSTER_ASSEMBLY,HSCB
# GOBP_2FE_2S_CLUSTER_ASSEMBLY,NFS1

# Specify the differential expression analysis you want to use for your gene
# list. This list will be filtered using adjusted p-value and fold change
# thresholds.
plot_id <- "f26be9c4-44b2-464e-a33c-66756f28b482"

# Define the adjusted p-value threshold to filter your gene list by.
adj_pvalue_threshold <- 0.01

# Define the fold change direction to filter your gene list by.
fold_change_direction <- "up" # up, down, or both

# Define the fold change threshold to filter your gene list by.
fold_change_threshold <- 1.2

# Define the minimum number of genes that a gene set must contain to be
# included in the analysis.
min_gene_set_size <- 5

# Define the maximum number of genes that a gene set must contain to be
# included in the analysis.
max_gene_set_size <- 1000

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Over-representation analysis was performed using the `clusterProfiler` R",
    " package. Passing genes were defined as those with an adjusted p-value of",
    " less than or equal to ", adj_pvalue_threshold, ", and a log fold-change",
    " greater than ", fold_change_threshold, "."
)

# Define plot colors for dots (colored by adjusted p-value).
low_color <- "blue"
high_color <- "red"

#------------------------------------------------------------------------------#
#----                          Helper function(s)                          ----#
#------------------------------------------------------------------------------#

get_gene_list <- function(
    df, adj_pvalue_threshold, fold_change_threshold, fold_change_direction) {
    if (fold_change_direction == "up") {
        filtered_df <- df %>%
            filter(
                Adj_P_Value <= adj_pvalue_threshold,
                Log2_Fold_Change > log2(fold_change_threshold)
            )
    } else if (fold_change_direction == "down") {
        filtered_df <- df %>%
            filter(
                Adj_P_Value <= adj_pvalue_threshold,
                Log2_Fold_Change < log2(fold_change_threshold)
            )
    } else if (fold_change_direction == "both") {
        filtered_df <- df %>%
            filter(
                Adj_P_Value <= adj_pvalue_threshold,
                Log2_Fold_Change < log2(fold_change_threshold) |
                    Log2_Fold_Change < -log2(fold_change_threshold)
            )
    }
    gene_list <- filtered_df$Gene_Symbol
    return(unique(gene_list))
}

format_ora_table <- function(ora_obj) {
    ora_table <- data.frame(ora_obj)
    ora_table <- ora_table[, -grep("Description", colnames(ora_table))]
    ora_table <- ora_table[order(ora_table$p.adjust, decreasing = FALSE), ]
    return(ora_table)
}

get_plot_data <- function(ora_table) {
    plot_data <- ora_table %>% slice_head(n = 20)
    plot_data <- plot_data %>%
        mutate(
            GeneRatio = sapply(GeneRatio, function(x) {
                parts <- strsplit(x, "/")[[1]]
                if (length(parts) == 2) {
                    as.numeric(parts[1]) / as.numeric(parts[2])
                } else {
                    NA # Handle cases where the format might be incorrect
                }
            })
        )
    # Reorder ID by GeneRatio
    plot_data$ID <- factor(
        plot_data$ID,
        levels = plot_data$ID[order(plot_data$GeneRatio)]
    )
    return(plot_data)
}

#------------------------------------------------------------------------------#
#----                              Main script                             ----#
#------------------------------------------------------------------------------#

# Load required libraries
library(pluto)
library(clusterProfiler)
library(tidyverse)

# Log into Pluto
pluto_login(api_token)

# Read in analysis results
deg_table <- pluto_read_results(
    experiment_id = experiment_id,
    plot_id = plot_id
)

# Get target gene sets
target_gene_sets <- read.csv(target_gene_set_file)

# Generate gene list from differential expression table
gene_list <- get_gene_list(
    deg_table, adj_pvalue_threshold,
    fold_change_threshold, fold_change_direction
)

# Run ORA
ora_obj <- enricher(
    gene_list,
    TERM2GENE = target_gene_sets,
    minGSSize = min_gene_set_size,
    maxGSSize = max_gene_set_size,
    qvalueCutoff = 1,
    pvalueCutoff = 1,
    universe = deg_table$Gene_Symbol,
    pAdjustMethod = "BH"
)

# Format ORA table
ora_table <- format_ora_table(ora_obj)

# Save ORA results
write.csv(ora_table, results_file_path, row.names = FALSE, na = "")

# Generate dot plot for top 20 pathways or gene sets
plot_data <- get_plot_data(ora_table)

plt <- ggplot(
    plot_data,
    aes(x = GeneRatio, y = ID, color = p.adjust, size = Count)
) +
    geom_point() +
    scale_color_gradient(low = low_color, high = high_color) +
    scale_size_continuous(range = c(2, 10)) +
    labs(
        x = "Gene Ratio",
        y = "Gene Set",
        color = "Adj. p-value",
        size = "Count"
    ) +
    theme_minimal() +
    theme(
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "in"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)
    )

ggsave(
    display_file_path,
    plot = plt,
    width = 12, height = 9, units = "in", dpi = 300
)

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)
