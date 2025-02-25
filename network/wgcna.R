################################################################################
#####             Weighted gene co-expression network analysis             #####
################################################################################

# Script to run weighted gene co-expression network analysis (WGCNA) on raw
# counts from a bulk RNA-seq experiment in Pluto Bio.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else,
# with the exception of updating the power_select parameter based on the output
# plots in Step 3.

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

# Define the phenotype column to use for the analysis.
# Note: this must match a column from your sample data!
phenotype_column <- "cell_line_modification"

# Define the trait to be tested for the analysis.
# Note: this must match a sample ID from your sample data!
tested_trait <- "R6 WT H3.3"

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

# Define a file path for the analysis .Rdat object.
analysis_object_file_path <- paste0(
  outdir, experiment_id, "_WGCNA_analysis.Rdat"
)

# Define a file path for the analysis plot (heatmap).
display_file_path <- paste0(
  outdir, experiment_id, "_WGCNA_heatmap.png"
)

# Include methods to describe your analysis.
plot_methods <- paste0(
  "Raw gene counts were filtered by removing outliers based on the function ",
  "goodSamplesGenes() from the WGCNA R package. Additionally, genes that did ",
  "not have at least ", min_counts, " counts in at least ", min_fraction * 100,
  "% of samples were removed from the analysis. Data was normalized using ",
  "variance stabilizing transformation (VST) from the DESeq2 R package. A ",
  "signed network was generated using a power threshold of ", power_select,
  " and using a Pearson correlation with a minimum module size of 30 genes. ",
  "All remaining parameters were used at default settings.",
  "All pairwise combinations of traits relative to all remaining samples, ",
  "excluding ", tested_trait, ", were binarized. Pearson correlations were ",
  "calculated for statistical significance within each comparison relative ",
  "to all remaining samples and plotted using the CorLevelPlot R package. ",
  "Asteriks indicate level of significance: * < 0.05, ** < 0.01, *** < 0.001, ",
  "while values indicate either a positive correlation (+) or an ",
  "anti-correlation (-) to the module eigengene."
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
#----        Step 3: Determine power value for network construction        ----#
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

# Plot graph
rsq / conn

#------------------------------------------------------------------------------#
#----                Step 4: Generate network based on power               ----#
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
  deepSplit = 2, # default
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
  file = analysis_object_file_path
)

#------------------------------------------------------------------------------#
#----               Step 5: Retrieve module eigengenes (MEs)               ----#
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
plotDendroAndColors(
  network$dendrograms[[1]],
  colors = cbind(network$unmergedColors, network$colors),
  c("modules before merge", "modules after merge"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

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
      outdir, experiment_id, "WGCNA_", me, "_eigengene_gene_members.txt",
      sep = ""
    )
  )
}

#------------------------------------------------------------------------------#
#----   Step 6: Associate module eigengene assignments to groups/traits    ----#
#------------------------------------------------------------------------------#
# Change reference level to tested_trait
all_levels <- unique(sample_data[[phenotype_column]])
levels_order <- c(tested_trait, setdiff(all_levels, tested_trait))
sample_data[[phenotype_column]] <- factor(
  sample_data[[phenotype_column]],
  levels = levels_order
)

# Binarize all associations
binarized_assocs <- binarizeCategoricalColumns(
  sample_data[[phenotype_column]],
  includePairwise = TRUE,
  includeLevelVsAll = TRUE,
  minCount = 1
)
colnames(binarized_assocs) <- gsub("^data\\.", "", colnames(binarized_assocs))

# Correlate each group versus with all remaining samples
binarized_assocs_vs_all <- binarized_assocs %>% select(ends_with(".vs.all"))

# Add sample names to rownames; should be the same order as sample_data
rownames(binarized_assocs_vs_all) <- rownames(sample_data)

# Get number of genes and samples for the all calculations
totalSamples <- nrow(binarized_assocs)
totalGenes <- ncol(wgcna_mat)

# Correlate the eigengenes to traits
# This step generates correlations and p-values, but note that it is not
# required for the visualization; this is just if you want the data calculated
# as a matrix. You could have also used a linear model or a logistic regression
# instead of a 0/1 correlation, but it does get you similar results.
eigengene_trait_assocs_ref_vs_all <- cor(
  module_eigengenes, binarized_assocs_vs_all,
  method = "pearson", use = "everything"
)
# The important thing to note here is that the association significance not
# necessarily the strength of the correlation:
eigengene_trait_assocs_pvals_ref_vs_all <- corPvalueStudent(
  eigengene_trait_assocs_ref_vs_all, totalSamples
)

#------------------------------------------------------------------------------#
#----                           Step 7: Visualize                          ----#
#------------------------------------------------------------------------------#
# Merge the eigengene info from step 5 with the binarized associations
vis_assoc_data_all <- merge(
  module_eigengenes, binarized_assocs_vs_all,
  by = "row.names"
)
rownames(vis_assoc_data_all) <- vis_assoc_data_all$Row.names

# Select only the columns you want to visualize
vis_assoc_data_all <- vis_assoc_data_all[, c(
  names(module_eigengenes), names(binarized_assocs_vs_all)
)]

# Generate a heatmap of the p-values and correlations
cor_plot <- CorLevelPlot(vis_assoc_data_all,
  x = names(binarized_assocs_vs_all), # column names of the traits
  y = names(module_eigengenes), # column names of the eigengene modules
  col = plot_colors,
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
png(display_file_path, width = 800, height = 600)
print(cor_plot)
dev.off()

#------------------------------------------------------------------------------#
#----                    Step 8: Push results to Pluto                     ----#
#------------------------------------------------------------------------------#
pluto_add_experiment_plot(
  experiment_id = experiment_id,
  display_file_path = display_file_path,
  analysis_name = analysis_name,
  plot_methods = plot_methods
)

#------------------------------------------------------------------------------#
#----                  Step 9: Find hub genes (optional)                   ----#
#------------------------------------------------------------------------------#
# Define column of binarized associations to explore
column_of_interest <- 1

# Calculate gene significance and associated p-values
gene_sig_corr <- cor(
  t(normalized_data), binarized_assocs_vs_all[[column_of_interest]],
  method = "pearson"
)

# Calculate the p-values
gene_sig_corr_pvals <- corPvalueStudent(gene_sig_corr, totalSamples)

# Get the top 25 hub genes that are significantly correlated to the association
hub_genes <- gene_sig_corr_pvals %>%
  as.data.frame() %>%
  arrange(V1) %>%
  head(25)
hub_genes$gene_symbol <- rownames(hub_genes)
colnames(hub_genes)[1] <- "pval"
hub_genes <- hub_genes[c("gene_symbol", "pval")]

write.csv(hub_genes,
  paste0(
    outdir,
    "hub_genes_for_",
    colnames(binarized_assocs_vs_all)[column_of_interest]
  ),
  row.names = FALSE,
  na = ""
)
