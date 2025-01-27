################################################################################
#####                     Gene set enrichment analysis                     #####
################################################################################

# Script to run gene set enrichment analysis (ORA) and create an enrichment plot
# or score bar plot of the top pathways or gene sets that are most significantly
# enriched in a given list of differentially expressed genes from Pluto Bio.

# This example script runs GSEA on a ranked list of genes. The genes should be
# ranked based on a metric such as fold change, p-value, or a combined ranking.
# To change the ranking method, modify the `rank_method` parameter under the
# user-defined parameters section.

# Additionally, this script highlights two visualizations you can generate from
# your GSEA results - an enrichment plot and a score bar plot. Note that we
# only define one display_file_path; this script is meant to be used such that
# you pick either the enrichment plot OR the score bar plot to generate and
# push back to Pluto (one plot per call to `pluto_add_experiment_plot()`).

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

# For more advanced plot customization, scroll down to the Main script and edit
# parameters for the  plot_ly function. For more information, see the usage
# docs: https://plotly.com/r/

################################################################################
#####                        User-defined parameters                       #####
################################################################################

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX185698"

# Speficy your experiment's organism.
# Options: human, mouse, rat, zebrafish, fruit fly
organism <- "mouse"

# Specify the collection shortname you would like to use from MSigDB.
# Note: to access collection shortnames programmatically, run
# unique(all_gene_sets$Collection_Shortname) after the prepare gene sets step.
collection_shortname <- "h"
collection_name <- "Hallmarks"

# Specify the differential expression analysis you want to use for your gene
# list. This list will be ranked according to rank_method.
plot_id <- "f393c752-1056-4772-b31c-76ac258a4006"

# Define the ranking method to use.
# Options: fold_change, p_value, sign_p_value, fold_change_p_value
rank_method <- "fold_change"

# Define the minimum number of genes that a gene set must contain to be
# included in the analysis.
min_gene_set_size <- 5

# Define the maximum number of genes that a gene set must contain to be
# included in the analysis.
max_gene_set_size <- 1000

# Define a name for the analysis in Pluto.
analysis_name <- "Gene set enrichment analysis"

# Define a file path for the analysis plot (volcano plot).
display_file_path <- paste0(experiment_id, "_", plot_id, "_gsea_plot.html")

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(experiment_id, "_", plot_id, "_gsea_results.csv")

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Gene set enrichment analysis (GSEA) was performed using the `fgsea` R ",
    "package and the `fgseaMultilevel()` function. The ", rank_method, " ",
    "from the differential expression comparison was used to rank genes. ",
    collection_name, " gene set collection from the Molecular Signatures ",
    "Database (MSigDB) was curated using the `msigdbr` R package. Prior to ",
    "running GSEA, the list of gene sets was filtered to include only gene ",
    "sets with between ", min_gene_set_size, " and ", max_gene_set_size, " ",
    "genes."
)

# Score bar plot parameters.
# Define plot colors for positively enriched and negatively enriched pathways.
pos_enriched_color <- "red"
neg_enriched_color <- "blue"
# Define the number of top pathways to display in the score bar plot.
top_pathways_num <- 10

################################################################################
#####                          Helper function(s)                          #####
################################################################################

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
    target_gene_sets <- split(
        x = target_gene_sets$gene_symbol,
        f = target_gene_sets$Gene_Set_Shortname
    )
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

get_ranked_gene_list <- function(deg_table, rank_method) {
    # Microarray data may have a suffix attached to duplicate gene names
    deg_table$Gene_Symbol <- paste0(deg_table$Gene_Symbol, "._x")
    deg_table$Gene_Symbol <- sapply(
        strsplit(deg_table$Gene_Symbol, split = "\\._"),
        function(x) {
            x[[1]]
        }
    )
    # Microarray and epigenetics data might have duplicate gene names
    # Keep the gene with the lowest p-value
    deg_table <- deg_table[!is.na(deg_table$P_Value), ]
    deg_table <- deg_table[order(deg_table$P_Value, decreasing = FALSE), ]
    deg_table <- deg_table[!(duplicated(deg_table$Gene_Symbol)), ]
    # If p-value is 0, set to min. non-zero value to avoid Inf values from log
    min_nonzero_p <- min(na.omit(deg_table$P_Value[deg_table$P_Value > 0]))
    deg_table[deg_table$P_Value == 0, "P_Value"] <- min_nonzero_p
    # Rank by selected method
    # Order by fold-change
    if (rank_method == "fold_change") {
        deg_table$rank_value <- deg_table$Log2_Fold_Change
        scoreType <- "std"
    }
    # Order by p-value
    if (rank_method == "p_value") {
        deg_table$rank_value <- -log10(deg_table$P_Value)
        scoreType <- "pos"
    }
    # Order by signed p-value
    if (rank_method == "sign_p_value") {
        deg_table$rank_value <- sign(deg_table$Log2_Fold_Change) *
            -log10(deg_table$P_Value)
        scoreType <- "std"
    }
    # Order by p-value * fold_change = the "lfcShrink" method
    if (rank_method == "fold_change_p_value") {
        deg_table$rank_value <- deg_table$Log2_Fold_Change *
            (-log10(deg_table$P_Value))
        scoreType <- "std"
    }
    # Generate ranked vector
    deg_table <- deg_table[order(deg_table$rank_value, decreasing = TRUE), ]
    deg_table$Rank <- 1:nrow(deg_table)
    ranks_vector <- as.numeric(deg_table$rank_value)
    names(ranks_vector) <- deg_table$Gene_Symbol
    return(list(
        ranked_gene_list = ranks_vector,
        ranked_deg_table = deg_table,
        score_type = scoreType
    ))
}

format_fgsea_table <- function(fgsea_table) {
    names(fgsea_table)[1] <- "Gene_Set_Shortname"
    # Order fgsea_table by significance
    fgsea_table <- fgsea_table[order(fgsea_table$padj, decreasing = FALSE), ]
    # Better formatting for the leading edge genes
    fgsea_table$leadingEdge <- as.character(fgsea_table$leadingEdge)
    fgsea_table$leadingEdge <- gsub(
        "[^A-Za-z0-9. ]", "", fgsea_table$leadingEdge
    )
    fgsea_table$leadingEdge <- gsub("^c", "", fgsea_table$leadingEdge)
    return(fgsea_table)
}

################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(fgsea)
library(msigdbr)
library(tidyverse)
library(plotly)
library(htmlwidgets)

#------------------------------------------------------------------------------#
#----                             Prepare data                             ----#
#------------------------------------------------------------------------------#
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

# Get ranked gene list from differential expression table
ranked_results <- get_ranked_gene_list(deg_table, rank_method)
ranked_gene_list <- ranked_results$ranked_gene_list
ranked_deg_table <- ranked_results$ranked_deg_table
score_type <- ranked_results$score_type

#------------------------------------------------------------------------------#
#----                    Gene set enrichment with fgsea                    ----#
#------------------------------------------------------------------------------#

# Run fgsea
fgsea_table <- fgseaMultilevel(
    pathways = target_gene_sets,
    stats = ranked_gene_list,
    minSize = min_gene_set_size,
    maxSize = max_gene_set_size,
    scoreType = score_type,
    eps = 0
)

# Format fgsea table
fgsea_table <- format_fgsea_table(fgsea_table)

# Save fgsea results
write.csv(fgsea_table, results_file_path, row.names = FALSE, na = "")

#------------------------------------------------------------------------------#
#----                       Create an enrichment plot                      ----#
#------------------------------------------------------------------------------#

# Note: to create an enrichment plot, you will need to know the name of the
# gene set that you want to plot. You can find the gene set names in the fgsea
# results table, under the Gene_Set_Shortname column.
gene_set <- "HALLMARK_MYOGENESIS"
enrichment_plot <- ggplotly(
    plotEnrichment(target_gene_sets[[gene_set]], ranked_gene_list)
)

# Save enrichment plot
saveWidget(
    enrichment_plot,
    file = display_file_path,
    title = "GSEA enrichment plot",
    selfcontained = TRUE
)

#------------------------------------------------------------------------------#
#----                        Create a score bar plot                       ----#
#------------------------------------------------------------------------------#

# Get top positively and negatively enriched pathways
top_pos <- fgsea_table %>%
    arrange(desc(NES)) %>%
    head(top_pathways_num)
top_neg <- fgsea_table %>%
    arrange(NES) %>%
    head(top_pathways_num)
top_pathways <- bind_rows(
    mutate(top_pos, direction = "Positively enriched"),
    mutate(top_neg, direction = "Negatively enriched")
)

# Modify the Gene_Set_Shortname column
top_pathways$Gene_Set_Shortname <- gsub(
    "^.*?_", "", top_pathways$Gene_Set_Shortname
)
top_pathways$Gene_Set_Shortname <- gsub(
    "_", " ", top_pathways$Gene_Set_Shortname
)

# Score bar plot for top pathways
score_bar_plot <- plot_ly(top_pathways,
    x = ~ reorder(Gene_Set_Shortname, (NES * -1)),
    y = ~NES,
    color = ~direction,
    colors = c(neg_enriched_color, pos_enriched_color),
    type = "bar"
) %>%
    layout(
        xaxis = list(
            title = "",
            tickangle = -90,
            showgrid = FALSE,
            zeroline = TRUE,
            showline = TRUE,
            linecolor = "black"
        ),
        yaxis = list(
            title = "Normalized enrichment score",
            showgrid = FALSE,
            zeroline = TRUE,
            showline = TRUE,
            linecolor = "black"
        ),
        barmode = "group",
        showlegend = TRUE,
        legend = list(
            traceorder = "reversed"
        )
    )

# Save score bar plot
saveWidget(
    score_bar_plot,
    file = display_file_path,
    title = "GSEA score bar plot",
    selfcontained = TRUE
)

#------------------------------------------------------------------------------#
#----                         Push results to Pluto                        ----#
#------------------------------------------------------------------------------#

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)
