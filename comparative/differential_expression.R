################################################################################
#####                   Differential expression analysis                   #####
################################################################################

# Script to compare genome-wide expression profiles of two groups and create a
# volcano plot visualization using data in a Pluto Bio experiment.

# Note that this script is designed for bulk RNA-seq data.

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
experiment_id <- "PLX185698"

# List the experimental sample IDs you want to use for the analysis.
# Note: these must match the sample IDs in your experiment sample data!
experimental_samples <- c("Fsk_24h_repA", "Fsk_24h_repB", "Fsk_24h_repC")
experiment_group_name <- c("Stimulated Fsk 24 hours")

# List the control sample IDs you want to use for the analysis.
# Note: these must match the sample IDs in your experiment sample data!
control_samples <- c("Ctrl_24h_repA", "Ctrl_24h_repB", "Ctrl_24h_repC")
control_group_name <- c("Control DMSO 24 hours")

# List the covariates you want to use for the analysis, NULL for none.
covariates <- NULL
# Covariates example: covariates <- c("sex")

# Assay data filtering options.
# Minimum number of reads a gene must have in a sample to be considered
read_count_threshold <- 3
# Minimum proportion of samples where the gene must meet `read_count_threshold`
sample_pct_threshold <- 0.2

# Define a file path for the analysis plot (volcano plot).
display_file_path <- paste0(
    experiment_id, "_deg_plot_",
    experiment_group_name, "_vs_", control_group_name,
    ".html"
)

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(
    experiment_id, "_deg_results_",
    experiment_group_name, "_vs_", control_group_name,
    ".csv"
)

# Define a name for the analysis in Pluto.
analysis_name <- "Differential expression analysis"

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Differential expression analysis was performed comparings the groups: ",
    experiment_group_name, " vs ", control_group_name, ". Genes were filtered ",
    "to include only genes with at least ", read_count_threshold, " reads ",
    "counted in at least ", sample_pct_threshold * 100, "% of samples in any ",
    "group. Differential expression analysis was then performed with the ",
    "DESeq2 R package, which tests for differential expression based on a ",
    "model using the negative binomial distribution."
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

validate_samples <- function(
    assay_data, sample_data, experimental_samples, control_samples) {
    assay_data_missing_ids <- setdiff(
        c(experimental_samples, control_samples), colnames(assay_data)
    )
    if (length(assay_data_missing_ids) > 0) {
        message(
            "Warning! The following samples are not in the assay data: ",
            paste(assay_data_missing_ids, collapse = ", ")
        )
    } else {
        message("All samples found in the assay data.")
    }
    sample_data_missing_ids <- setdiff(
        c(experimental_samples, control_samples), sample_data$sample_id
    )
    if (length(sample_data_missing_ids) > 0) {
        message(
            "Warning! The following samples are not in the sample data: ",
            paste(sample_data_missing_ids, collapse = ", ")
        )
    } else {
        message("All samples found in the sample data.")
    }
}

format_assay_data <- function(assay_data) {
    # Ensure unique gene names in assay data
    assay_data <- assay_data %>%
        mutate(gene_symbol = make.unique(gene_symbol, sep = "_")) %>%
        column_to_rownames("gene_symbol")
    # Check for non-integer columns and convert if needed
    non_integer_columns <- assay_data %>%
        summarise_all(~ all(. %% 1 == 0)) %>%
        gather() %>%
        filter(value == FALSE) %>%
        pull(key)
    if (length(non_integer_columns) > 0) {
        message(
            "Non-integer columns detected: ",
            paste(non_integer_columns, collapse = ", ")
        )
        message("Double-check those columns to make sure they are correct!")
        assay_data <- assay_data %>%
            mutate(across(everything(), as.integer))
    }
    message(
        "Assay data contains ", nrow(assay_data),
        " genes and ", ncol(assay_data), " samples total"
    )
    return(assay_data)
}

subset_assay_data <- function(
    assay_data, experimental_samples, control_samples) {
    assay_data <- assay_data %>%
        select(
            intersect(experimental_samples, names(assay_data)),
            intersect(control_samples, names(assay_data))
        )
    message(
        "Assay data dimensions: ", nrow(assay_data),
        " genes x ", ncol(assay_data), " samples"
    )
    return(assay_data)
}

filter_assay_data <- function(
    assay_data, experimental_samples, control_samples,
    sample_pct_threshold, read_count_threshold) {
    sample_threshold <- length(
        c(experimental_samples, control_samples)
    ) * sample_pct_threshold
    genes_to_keep <- apply(assay_data, 1, function(x) {
        sum(x >= read_count_threshold) >= sample_threshold
    })
    assay_data <- assay_data[genes_to_keep, ]
    message(nrow(assay_data), " genes left after filtering")
    return(assay_data)
}

create_deseq2_sample_data <- function(
    sample_data, experimental_samples, control_samples, covariates) {
    # Match experimental and control samples with experiment sample data
    sample_data <- sample_data[match(
        c(experimental_samples, control_samples),
        sample_data$sample_id
    ), ]
    # Confirm that all covariates are column names in the sample data
    if (all(
        !is.null(covariates) & !all(covariates %in% colnames(sample_data))
    )) {
        stop("The following covariates are not in sample data input: ",
            covariates[!(covariates %in% colnames(sample_data))],
            call. = FALSE
        )
    }
    # Create sample data for DESeq2 and combine with covariates
    deseq2_sample_data <- cbind(
        data.frame(
            sample_id = c(experimental_samples, control_samples),
            condition = c(
                rep("experiment", length(experimental_samples)),
                rep("control", length(control_samples))
            )
        ),
        sample_data[, covariates]
    )
    if (ncol(deseq2_sample_data) > 2) {
        colnames(deseq2_sample_data)[3:ncol(deseq2_sample_data)] <- covariates
    }
    return(deseq2_sample_data)
}

run_deseq2 <- function(sample_data, assay_data) {
    message("Running differential expression with DESeq2...")
    message(paste0(
        "Comparing experimental group (n = ", length(experimental_samples),
        ") vs control group (n = ", length(control_samples), ")"
    ))
    if (ncol(sample_data) > 2) {
        design <- paste0(
            "~ ", paste0(
                colnames(sample_data)[3:ncol(sample_data)],
                collapse = " + "
            ),
            " + condition"
        )
    } else {
        design <- "~ condition"
    }
    dds <- DESeq(DESeqDataSetFromMatrix(
        countData = assay_data,
        colData = sample_data,
        design = as.formula(design)
    ))
    res <- results(dds,
        contrast = c("condition", "experiment", "control"),
        pAdjustMethod = "fdr",
        cooksCutoff = FALSE,
        independentFiltering = FALSE
    )
    return(res)
}

get_deseq2_table <- function(res) {
    deg_table <- data.frame(
        res[, c("baseMean", "log2FoldChange", "pvalue", "padj")],
        stringsAsFactors = FALSE
    )
    deg_table <- cbind(
        Gene_Symbol = row.names(deg_table), deg_table, stringsAsFactors = FALSE
    )
    deg_table <- deg_table[order(deg_table$padj, decreasing = FALSE), ]
    # Clean up column names
    names(deg_table) <- c(
        "Gene_Symbol", "Average_Expression", "Log2_Fold_Change",
        "P_Value", "Adj_P_Value"
    )
    return(deg_table)
}

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
library(DESeq2)
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

# Read in experiment sample data
sample_data <- pluto_read_sample_data(experiment_id)

# Validate experimental and control samples
validate_samples(assay_data, sample_data, experimental_samples, control_samples)

#------------------------------------------------------------------------------#
#----                  Differential expression with DESeq2                 ----#
#------------------------------------------------------------------------------#
# Format experiment assay data for DESeq2
deseq2_assay_data <- format_assay_data(assay_data)

# Subset assay data to just the control and experimental samples
deseq2_assay_data <- subset_assay_data(
    deseq2_assay_data, experimental_samples, control_samples
)

# Filter assay data to genes with read_count_threshold reads in
# n*sample_pct_threshold samples
deseq2_assay_data <- filter_assay_data(
    deseq2_assay_data, experimental_samples, control_samples,
    sample_pct_threshold, read_count_threshold
)

# Create sample data for DESeq2
deseq2_sample_data <- create_deseq2_sample_data(
    sample_data, experimental_samples, control_samples, covariates
)

# Run DESeq2
deseq2_res <- run_deseq2(deseq2_sample_data, deseq2_assay_data)

# Get results table
deg_table <- get_deseq2_table(deseq2_res)
write.csv(deg_table, results_file_path, row.names = FALSE)

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
    title = paste0(experiment_group_name, " vs ", control_group_name),
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
