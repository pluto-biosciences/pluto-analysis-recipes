################################################################################
#####          Weighted gene co-expression network analysis (all)          #####
################################################################################

# Script to run weighted gene co-expression network analysis (WGCNA) on raw
# counts from a bulk RNA-seq experiment in Pluto Bio.

# In this version of WGCNA, all trait associations are binarized. This is the
# standard / traditional way to perform WGCNA.
# i.e., compare all trait levels against each other

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else,
# with the exception of updating the power_select parameter based on the output
# plots in Step 2.

# For more advanced plot customization, scroll down to the Main script and edit
# parameters for the  CorLevelPlot function. For more information, see the usage
# docs: https://github.com/kevinblighe/CorLevelPlot

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

# Define the trait to be tested for the analysis.
# Note: this must match a sample ID from your sample data!
tested_trait <- "cell_line_modification"

# Assay data filtering options.
# Minimum number of counts per gene per sample to be included in the analysis
min_counts <- 3
# Minimum fraction of samples that need to have min_counts to retain a gene
min_fraction <- 0.20

# WGCNA parameters.
# Soft-thresholding power for network construction
# Note: this value should be updated based on power analysis results!
power_select <- 12
# Network type for WGCNA analysis ("signed", "unsigned", or "signed hybrid")
# "unsigned" means the direction of correlation does not matter; genes are
# connected in a network whether it is a positive or negative correlation.
# "signed" means the direction of correlation does matter; whether a gene is
# positively or negatively correlated is taken into account when deciding on a
# connection between genes.
network_type <- "signed"
# Number of threads to use for parallel processing
threads <- 10

# Define a name for the analysis in Pluto.
analysis_name <- "WGCNA"

# Define an outdir for WGCNA results.
outdir <- paste0(
    experiment_id, "_WGCNA_analysis/"
)
# Create the directory if it doesn't exist
if (!dir.exists(outdir)) {
    dir.create(outdir)
}

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Raw gene counts were filtered by removing outliers based on the function ",
    "goodSamplesGenes() from the WGCNA R package. Additionally, genes that ",
    "did not have at least ", min_counts, " counts in at least ",
    min_fraction * 100,
    "% of samples were removed from the analysis. Data was normalized using ",
    "variance stabilizing transformation (VST) from the DESeq2 R package. A ",
    "signed network was generated using a power threshold of ", power_select,
    " and using a Pearson correlation with a minimum module size of 30 genes. ",
    "All remaining parameters were used at default settings.",
    "All pairwise combinations of ", tested_trait, " relative to all ",
    "remaining samples were binarized. Pearson correlations were ",
    "calculated for statistical significance within each comparison relative ",
    "to all remaining samples and plotted using the CorLevelPlot R package. ",
    "Asteriks indicate level of significance: * < 0.05, ** < 0.01, ",
    "*** < 0.001, while values indicate either a positive correlation (+) or ",
    "an anti-correlation (-) to the module eigengene."
)

# Define plot colors for heatmap.
plot_colors <- c("blue1", "skyblue", "white", "pink", "red")

################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(WGCNA)
library(DESeq2)
library(tidyverse)
library(patchwork)
library(assertthat)
library(CorLevelPlot)
allowWGCNAThreads() # Allow multi-threading for running WGCNA functions

#------------------------------------------------------------------------------#
#----                         Step 0: Prepare data                         ----#
#------------------------------------------------------------------------------#
# Log into Pluto
pluto_login(api_token)

# Read in experiment assay data
assay_data <- pluto_read_assay_data(experiment_id)
rownames(assay_data) <- assay_data$gene_symbol
assay_data <- assay_data[, -1]

# Read in experiment sample data
sample_data <- pluto_read_sample_data(experiment_id)
rownames(sample_data) <- sample_data$sample_id
sample_data <- sample_data[, -1]
sample_data <- sample_data[colnames(assay_data), ]

# Confirm variables
assert_that(is.count(min_counts))
assert_that(((min_fraction >= 0) & (min_fraction <= 1)))
assert_that(is.count(power_select))
assert_that(network_type %in% c("signed", "unsigned", "signed hybrid"))
assert_that(tested_trait %in% colnames(sample_data))

timestamp <- str_replace_all(
    string = Sys.time(), pattern = " ", replacement = "_"
)

#------------------------------------------------------------------------------#
#----                  Step 1: Filter and normalize data                   ----#
#------------------------------------------------------------------------------#
# Exclude any genes or samples that are outliers
check_for_outliers <- WGCNA::goodSamplesGenes(t(assay_data))
table(check_for_outliers$goodGenes)
assay_data <- assay_data[check_for_outliers$goodGenes == TRUE, ]

# Filter genes and only those that have at least min_counts counts in at least
# min_fraction of the samples
filtered_assay_data <- assay_data[
    names(which(rowSums(
        assay_data[, ] >= min_counts
    ) >= round(ncol(assay_data) * min_fraction))),
]

# Normalize data using DESeq2 variance stabilizing transform
# Note: this is not running DE analysis, so no model is necessary. We are just
# getting normalized counts (via variable stabilizing transformation).
de_obj <- DESeqDataSetFromMatrix(round(filtered_assay_data),
    sample_data,
    design = ~1
)
de_norm_obj <- vst(de_obj)
normalized_data <- assay(de_norm_obj)

# Transpose data for WGCNA
wgcna_mat <- t(normalized_data)

#------------------------------------------------------------------------------#
#----        Step 2: Determine power value for network construction        ----#
#------------------------------------------------------------------------------#
# In the pickSoftThreshold() function, the powerVector parameter is the default
# provided by the WGCNA package. You can provide a different list/range of
# powers to test in order to find the best power for your analysis.
select_soft_thresh <- pickSoftThreshold(
    data = wgcna_mat,
    dataIsExpr = TRUE,
    powerVector = c(seq(1, 10, by = 1), seq(12, 50, by = 2)),
    RsquaredCut = 0.85, # default
    nBreaks = 10, # default
    networkType = network_type,
    corFnc = cor, # default
    corOptions = list(use = "p"), # default
    verbose = 5
)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# !!!    You will need to update power_select based on the graphs below    !!! #
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! #
# Plot the data to find what power threshold to use; similar to an elbow plot,
# pick a power near the elbow curve. Aim to pick a threshold that maximizes
# Rsquared while minimizing mean K connectivity.

# Power for r-squared (pick something above the line, but don't overfit)
rsq <- ggplot(
    select_soft_thresh$fitIndices, aes(x = Power, y = SFT.R.sq, label = Power)
) +
    geom_point() +
    geom_text(nudge_y = 0.05) +
    geom_hline(yintercept = 0.8, color = "orange", linetype = "dashed") +
    labs(
        x = "Power",
        y = paste(
            "Scale free toplogy model fit,", network_type, "R-squared",
            sep = " "
        )
    ) +
    ggtitle("Power determination for network construction") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

# Power for mean connectivity (smaller is better, but don't overfit)
conn <- ggplot(
    select_soft_thresh$fitIndices, aes(x = Power, y = mean.k., label = Power)
) +
    geom_point() +
    geom_text(nudge_x = 0.7) +
    labs(x = "Power", y = "mean connectivity") +
    theme_minimal() +
    ggtitle("")

# Plot
rsq / conn

# Save plot
power_plot <- rsq / conn
ggsave(filename = paste(
    outdir, "power_determination_for_network_", timestamp, ".png",
    sep = ""
), device = "png", plot = power_plot, width = 10, height = 6)

#------------------------------------------------------------------------------#
#----                Step 3: Generate network based on power               ----#
#------------------------------------------------------------------------------#
# Save current cor() namespace to temp var since will use WGCNA cor()
temp_cor <- cor
# Force it to use WGCNA cor() (fix a namespace conflict issue)
cor <- WGCNA::cor

# Call the network topology analysis function
network <- blockwiseModules(
    datExpr = wgcna_mat, power = power_select, networkType = network_type,
    corType = "pearson", # default
    TOMType = network_type,
    mergeCutHeight = 0.25, # threshold to use to merge similar modules; default
    numericLabels = FALSE, # modules will be color names instead of numbers
    deepSplit = 1, # default
    detectCutHeight = 0.995, # default
    minModuleSize = 30, # default: min(20, ncol(datExpr)/2 )
    maxBlockSize = 22000, # default: 5000, based on available computer memory
    randomSeed = 42,
    nThreads = threads,
    useCorOptionsThroughout = TRUE,
    verbose = 5
)
cor <- temp_cor # Return cor() to original namespace

# Save WGCNA analysis to .Rdat object
save(
    list = c(
        "select_soft_thresh",
        "wgcna_mat",
        "filtered_assay_data",
        "normalized_data",
        "sample_data",
        "network"
    ),
    file = paste(
        outdir, experiment_id, "_WGCNA_analysis_", timestamp, ".Rdat",
        sep = ""
    )
)

#------------------------------------------------------------------------------#
#----               Step 4: Retrieve module eigengenes (MEs)               ----#
#------------------------------------------------------------------------------#
# Retrieve module eigengenes (MEs)
module_eigengenes <- network$MEs

# Number of genes in each MEs
table(network$colors)

# Plot dendrograms and module colors before and after merging
# Merging means some modules were very similar and if they hit the merge
# threshold listed above then they were merged into one module. Use the merged
# modules for further analysis. Grey modules are genes that don't fit in any of
# the modules, so they are assigned to the grey category.
png(filename = paste(
    outdir, "dendrogram_before_after_module_merge_", timestamp, ".png",
    sep = ""
), width = 900, height = 700)
dendro_merge <- plotDendroAndColors(
    network$dendrograms[[1]],
    colors = cbind(network$unmergedColors, network$colors),
    c("modules before merge", "modules after merge"),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
)
dev.off()

par(mar = c(5.1, 4.1, 4.1, 2.1))
png(filename = paste(
    outdir, "correlation_of_eigengene_module_similarity_", timestamp, ".png",
    sep = ""
), width = 1000, height = 700)
if (network_type == "signed") {
    WGCNA::plotEigengeneNetworks(
        multiME = network$MEs,
        setLabels = colnames(module_eigengenes),
        excludeGrey = TRUE,
        signed = TRUE,
        plotAdjacency = FALSE,
        plotDendrograms = TRUE,
        colorLabels = TRUE,
        plotHeatmaps = TRUE,
        setMargins = FALSE
    )
} else {
    WGCNA::plotEigengeneNetworks(
        multiME = network$MEs,
        setLabels = colnames(module_eigengenes),
        excludeGrey = TRUE,
        signed = FALSE,
        plotAdjacency = FALSE,
        plotDendrograms = TRUE,
        colorLabels = TRUE,
        plotHeatmaps = TRUE,
        setMargins = FALSE
    )
}
dev.off()

par(mar = c(5.1, 4.1, 4.1, 2.1))
png(filename = paste(
    outdir, "relationship_of_eigengene_modules_and_samples_", timestamp, ".png",
    sep = ""
), width = 1200, height = 1200)
WGCNA::plotMEpairs(
    datME = module_eigengenes[
        colnames(module_eigengenes)[
            which(!colnames(module_eigengenes) %in% "MEgrey")
        ]
    ], clusterMEs = TRUE
)
dev.off()

# Get the genes associated with the module colors based on the heatmap
genes_in_modules <- as.data.frame(network$colors)
names(genes_in_modules)[1] <- "MEcolor"
genes_in_modules$genes <- rownames(genes_in_modules)
# Split genes into MEs for each file
me_genes <- split(genes_in_modules$genes, genes_in_modules$MEcolor)
for (me in names(me_genes)) {
    writeLines(
        text = me_genes[[me]],
        con = paste(
            outdir, "WGCNA_", me, "_eigengene_gene_members_", timestamp, ".txt",
            sep = ""
        )
    )
}

#------------------------------------------------------------------------------#
#----   Step 5: Associate module eigengene assignments to groups/traits    ----#
#------------------------------------------------------------------------------#
# Identifies modules that are associated to traits
# Categorize traits by binarization to find associations
# Creates a separate column for each level

# Convert to factor/categorical
sample_data[, tested_trait] <- factor(sample_data[, tested_trait])

# Binarize all associations
binarized_assocs <- binarizeCategoricalColumns(
    sample_data[, tested_trait],
    levelSep.vsAll = " vs ", includeLevelVsAll = TRUE, minCount = 1
)
colnames(binarized_assocs) <- gsub("^data\\.", "", colnames(binarized_assocs))

# Add sample names to rownames; should be the same order as sample_data
rownames(binarized_assocs) <- rownames(sample_data)

# Correlate the eigengenes to traits
# This step generates correlations and p-values, but note that it is not
# required for the visualization; this is just if you want the data calculated
# as a matrix. You could have also used a linear model or a logistic regression
# instead of a 0/1 correlation, but it does get you similar results.
totalSamplesInCorr <- nrow(binarized_assocs)
eigengene_trait_assocs <- cor(
    module_eigengenes, binarized_assocs,
    method = "pearson", use = "all"
)
eigengene_trait_assocs_pvals <- corPvalueStudent(
    eigengene_trait_assocs, totalSamplesInCorr
)

#------------------------------------------------------------------------------#
#----                           Step 6: Visualize                          ----#
#------------------------------------------------------------------------------#
# Merge the eigengene info from step 4 with the binarized associations
vis_assoc_data_all <- merge(
    module_eigengenes, binarized_assocs,
    by = "row.names"
)
rownames(vis_assoc_data_all) <- vis_assoc_data_all$Row.names

# Select only the columns you want to visualize
vis_assoc_data_all <- vis_assoc_data_all[, c(
    names(module_eigengenes), names(binarized_assocs_vs_all)
)]

# Generate a heatmap of the p-values and correlations
cor_plot <- CorLevelPlot(vis_assoc_data_all,
    x = names(binarized_assocs), # column names of the traits
    y = names(module_eigengenes), # column names of the eigengene modules
    col = plot_colors,
    corFUN = "pearson",
    main = "ME associations with binarized phenotype trait comparisons",
    cexMain = 2.0,
    titleY = "Module Eigengene (ME)",
    rotTitleY = 90,
    cexTitleY = 1.5,
    fontLabY = 1,
    rotLabY = 30,
    titleX = "Binarized traits",
    cexTitleX = 1.5,
    fontLabX = 1,
    rotLabX = 45
)

# Save heatmap
display_file_path <- paste(
    outdir, "correlation_of_eigengene_modules_with_binarized_traits_",
    timestamp, ".png",
    sep = ""
)
png(display_file_path, width = 800, height = 600)
print(cor_plot)
dev.off()

#------------------------------------------------------------------------------#
#----                    Step 7: Push results to Pluto                     ----#
#------------------------------------------------------------------------------#
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)

#------------------------------------------------------------------------------#
#----                  Step 8: Find hub genes (optional)                   ----#
#------------------------------------------------------------------------------#
# Find hub genes in the modules (genes with high module membership)
module_membership_strength <- cor(
    module_eigengenes, t(normalized_data),
    method = "pearson"
)
totalSamples <- nrow(t(normalized_data))
module_membership_strength_pvals <- corPvalueStudent(
    module_membership_strength, totalSamples
)

write.csv(x = t(module_membership_strength), file = paste(
    outdir, network_type, "_pearson_correlation_weights_for_each_gene_contributing_to_each_eigengene_module_", timestamp, ".csv",
    sep = ""
), row.names = TRUE, col.names = TRUE)
write.csv(x = t(module_membership_strength_pvals), file = paste(
    outdir, network_type, "_pearson_correlation_pvalues_for_each_gene_contributing_to_each_eigengene_module_", timestamp, ".csv",
    sep = ""
), row.names = TRUE, col.names = TRUE)

# Get correlation and p-values for comparisons associated with genes
# In particular, when using the "all" association, this can be used as sample
# size does not change.
gene_sig_corr <- cor(t(normalized_data), binarized_assocs, method = "pearson")
gene_sig_corr_pvals <- corPvalueStudent(gene_sig_corr, totalSamples)

top_genes_to_plot <- c()
for (comp in colnames(binarized_assocs)) {
    print(paste("Running", comp))
    correlations <- gene_sig_corr[, comp, drop = FALSE]
    pvalues <- gene_sig_corr_pvals[, comp, drop = FALSE]
    colnames(correlations) <- "pearson_correlation_coeff"
    colnames(pvalues) <- "pvalue"
    assocs <- merge(correlations, pvalues, by = "row.names")
    colnames(assocs)[1] <- "gene"
    top_genes_in_me <- assocs[ # Select top genes by lowest pvalues
        order(assocs$pvalue, decreasing = FALSE)[1:25], "gene"
    ]
    names(top_genes_in_me) <- rep(x = comp, length(top_genes_in_me))
    top_genes_to_plot <- c(top_genes_to_plot, top_genes_in_me)
    write.csv(x = assocs, file = paste(
        outdir, comp, "_gene_correlation_pval_matrix_", timestamp, ".csv",
        sep = ""
    ), row.names = FALSE, col.names = TRUE)
}

# Plot heatmap of network connectivity of genes
WGCNA::plotNetworkHeatmap(
    datExpr = wgcna_mat,
    plotGenes = top_genes_to_plot,
    useTOM = FALSE,
    power = power_select,
    networkType = network_type
)
