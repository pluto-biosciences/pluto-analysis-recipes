################################################################################
#####                    Overlap analysis - peaks calls                    #####
################################################################################

# Script to create a Venn diagram of overlapping peaks from multiple
# independent peaks calls for each category in Pluto Bio.

# This example script creates a Venn diagram for overlapping peaks across
# categories, which can be either individual samples or consensus peaks.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

# For more advanced plot customization, scroll down to the Main script and edit
# parameters for the findOverlapsOfPeaks() and makeVennDiagram() functions.
# For more information, see the usage docs:
# https://bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/
# inst/doc/ChIPpeakAnno.html#42_Assess_the_concordance_of_peak_replicates

################################################################################
#####                        User-defined parameters                       #####
################################################################################

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX075908"

# Define a file path for the analysis plot (Venn diagram).
display_file_path <- paste0(experiment_id, "_venn_diagram.png")

# Define a name for the analysis in Pluto.
analysis_name <- "Overlap analysis: Peaks"

# Include methods to describe your analysis.
plot_methods <- paste(
  "Overlapping peaks across samples,",
  "plotted using the `R` `ChIPpeakAnno` package."
)

# List the .bed files to use for overlap analysis (at least 2, less than 5).
# You can compare consensus peak .bed files across experiments, or you can
# compare individual samples within/across experiments. Reach out to Pluto
# support for help accessing your .bed files.
bed_files <- c(
  "AP-1.consensus_peaks.annotatePeaks.bed",
  "NF-kB.consensus_peaks.annotatePeaks.bed",
  "Sp1.consensus_peaks.annotatePeaks.bed"
)

# List category names that correspond to each of the above .bed files.
# Note: order matters! The order of category_names should match the order of
# bed_files.
category_names <- c("AP-1", "NF-kB", "Sp1")

# Define plot colors (one color per .bed file / category name).
# Note: the length of plot_colors must match the length of bed_files and
# the length of category_names.
plot_colors <- c("cornflowerblue", "orchid1", "pink")

# Define peaks overlapping criteria between samples.
# Parameter definitions:
#     * maxgap: default is -1L which means a peak will be overlapping if one
#       range has its start or end strictly inside the other; otherwise, setting
#       this parameter to a value such as 20 means that two peaks will be
#       considered overlapping if the number of base pairs that separates them
#       is less than 20 base pairs apart.
#     * minoverlap: default is 0L. This is a value between 0.0-1.0. A peak is
#       considered overlapping if the shared fraction is at least the values
#       indicated. For example, if this was set to 0.10, then peaks would only
#       be considered overlapping if at least 10% of the base pairs are shared
#       between them.
#     * connectedPeaks: options are keepAll or keepFirstListConsistent; this is
#       because multiple peaks in one category might map to one broad peak in
#       another category. Setting "keepAll" will tell you how many peaks overlap
#       from each category, while setting "keepFirstListConsistent" will keep
#       the number in the Venn diagram based on the first category you listed
#       above (Note: you won't get the peak number break down per category).
maxgap <- -1L
minoverlap <- 0L # must be between 0.0 and 1.0
connectedPeaks <- "keepAll" # options: "keepAll" or "keepFirstListConsistent"

################################################################################
#####                              Main script                             #####
################################################################################

# Load required libraries
library(pluto)
library(tidyverse)
library(ChIPpeakAnno)
library(grid)

# Log into Pluto
pluto_login(api_token)

# Convert .bed files to GRanges objects
peaks_regions_per_category <- lapply(bed_files, function(bedFile) {
  toGRanges(bedFile, format = "BED", header = FALSE)
})

# Name the Granges objects to match sample names
names(peaks_regions_per_category) <- category_names

# Get peak overlaps
# Note: the connectedPeaks in findOverlapsOfPeaks() is not the same parameter as
# connectedPeaks in makeVennDiagram()
overlapping_peaks <- findOverlapsOfPeaks(
  peaks_regions_per_category,
  maxgap = maxgap,
  minoverlap = minoverlap,
  ignore.strand = TRUE,
  connectedPeaks = "keepAll"
)

# Plot Venn diagram for all peaks in each category
png(
  display_file_path,
  width = 1200,
  height = 1200,
  res = 300
)

venn_plot <- makeVennDiagram(
  overlapping_peaks,
  fill = plot_colors,
  connectedPeaks = connectedPeaks,
  plot = TRUE,
  lwd = 2,
  lty = "blank",
  cex = 0.75, # change font size of overlap values
  fontfamily = "sans", # change font family of overlap values
  cat.cex = 1.2, # change font size of category names
  cat.fontface = "bold", # change font style of category names
  cat.fontfamily = "sans", # change font family of category names
  cat.default.pos = "outer"
)

dev.off()

# Push results back to Pluto
pluto_add_experiment_plot(
  experiment_id = experiment_id,
  display_file_path = display_file_path,
  analysis_name = analysis_name,
  plot_methods = plot_methods
)

# Note: the overlapping_peaks object contains a lot of information that might
# be interesting depending on your analysis questions:
#     * venn_cnt: an object of VennCounts
#     * peaklist: a list consists of all overlapping peaks or unique peaks
#     * uniquePeaks: an object of GRanges consists of all unique peaks
#     * mergedPeaks: an object of GRanges consists of all merged overlapping
#       peaks
#     * peaksInMergedPeaks: an object of GRanges consists of all peaks in each
#     * category involved in the overlapping peaks
#     * overlappingPeaks: a list of data frame consists of the annotation of all
#       the overlapped peaks
#     * all.peaks: a list of GRanges object which contain the input peaks with
#       formated rownames.
