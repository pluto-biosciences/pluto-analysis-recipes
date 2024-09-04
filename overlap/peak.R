################################################################################
#####           Overlap analysis - independent sample peaks calls          #####
################################################################################

# Script to create a Venn diagram of overlapping peaks from multiple
# independent peaks calls for each sample in Pluto Bio.

# This example script creates a Venn diagram for overlapping peaks across 
# individual samples.

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
experiment_id <- "PLX278874"

# Define a file path for the analysis plot (Venn diagram).
display_file_path <- paste0(experiment_id, "_venn_diagram.png")

# Define a name for the analysis in Pluto.
analysis_name <- "Overlap analysis: Peaks"

# Include methods to describe your analysis.
plot_methods <- paste(
  "Overlapping peaks across samples,",
  "plotted using the `R` `ChIPpeakAnno` package."
)

# List the samples you want to include in your peak analysis.
# The sample ids must match the sample names listed your Pluto sample data.
sample_ids <- c("sample1", "sample2", "sample3")

# List sample names that correspond to each of the above analyses.  This will be
# what is displayed on the graph and does not need to match the Pluto sample
# data sample ids. Note: order matters! The order of sample_names should match
# the order of sample_ids.
sample_names <- c("H3K27ac", "H3K27me3", "TF")

# Define venn circle fill and border outline colors (one color per sample).
# Note: the length of venn_circle_colors and venn_circle_border_colors must
# match the length of samples_ids and sample_names.
venn_circle_colors <- c("cornflowerblue", "orchid1", "pink")
venn_circle_border_colors <- c("blue4", "darkviolet", "deeppink")

# Define peaks overlapping criteria between samples.
# Parameter definitions:
#     * maxgap: default is -1L which means a peak will be overlapping if one
#     range has its start or end strictly inside the other; otherwise, setting
#     this parameter to a value such as 20 means that two peaks will be
#     considered overlapping if the number of base pairs that separates them is
#     less than 20 base pairs apart.
#     * minoverlap:  default is 0L.  This is a value between 0.0-1.0.  A peak is
#     considered overlapping if the shared fraction is at least the values 
#     indicated.  For example, if this was set to 0.10, then peaks would only
#     be considered overlapping if at least 10% of the base pairs are shared
#     between them.
#     * connectedPeaks: options are keepAll or keepFirstListConsistent; this is 
#     because multiple peaks in one sample might map to one broad peak in
#     another sample.  If you set "keepAll" then it will tell you how many 
#     peaks overlap from each sample, while setting "keepFirstListConsistent"
#     means to keep the number in the venn diagram based on the first sample
#     you listed above (Note: you won't get the peak number break down per 
#     sample).
maxgap = -1L
minoverlap = 0L
connectedPeaks = "keepAll"

################################################################################
#####                          Helper function(s)                          #####
################################################################################

# Confirm all input values are valid
validate_input_parameters <- function(){
  if (length(sample_ids) != length(sample_names)){
    stop("The total number of samples listed (", length(sample_ids), ") does
         not match the total number of sample names (", length(sample_names),
         ").")
  }
  if (length(sample_ids) !=  length(venn_circle_colors)){
    stop("The total number of samples listed(", length(sample_ids), ") does not
         match the total number of venn_circle_colors listed (", 
         length(venn_circle_colors), ").")
  }
  if (length(sample_ids) !=  length(venn_circle_border_colors)){
    stop("The total number of samples listed(", length(sample_ids), ") does not
         match the total number of venn_circle_border_colors listed (", 
         length(venn_circle_border_colors), ").")
  }
  if (minoverlap > 1 | minoverlap < 0){
    stop("The minoverlap parameter must be between 0.0 and 1.0")
  }
  if(connectedPeaks != "keepAll" & connectedPeaks != "keepFirstListConsistent"){
    stop("The connectedPeaks parameters can only be keepAll or 
         keepFirstListConsistent.")
  }
}

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

# Confirm all user defined inputs are valid
validate_input_parameters()

# Convert per-sample peak bed files to Granges objects
peaks_regions_per_sample <- lapply(results_list, function(bedFile) {
  toGRanges(bedFile, format = "BED", header = FALSE)
  }
)

# Name the Granges objects to match sample names
names(peaks_regions_per_sample) <- sample_names

# Get peak overlaps
# Note: the connectedPeaks in findOverlapsOfPeaks() is not the same parameter as
# connectedPeaks in makeVennDiagram()
overlapping_peaks <- findOverlapsOfPeaks(
  peaks_regions_per_sample, maxgap = maxgap, minoverlap = minoverlap, 
  ignore.strand = TRUE, connectedPeaks = "keepAll"
)

# Plot Venn diagram for all peaks in each sample
venn_plot <- makeVennDiagram(
  overlapping_peaks, fill = venn_circle_colors, col = venn_circle_border_colors,
  connectedPeaks = connectedPeaks, plot = TRUE
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