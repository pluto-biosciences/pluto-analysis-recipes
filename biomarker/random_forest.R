################################################################################
#####                 Random forest for biomarker selection                #####
################################################################################

# This script runs a Random Forest model to perform biomarker selection based
# on assay data and sample data from a Pluto Bio experiment, using a selected
# categorical variable from the sample data and specifying two groups for
# comparison. The analysis will then identify key biomarkers distinguishing
# these groups.

# This script will output a decision tree plot based on a CART model using the
# top features identified by the Random Forest model as being important,
# focusing on the most important biomarkers identified across the Random Forest.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and variable. You should not need to change anything else.

# For more advanced configuration options and fine-tuning of the Random Forest
# model, scroll down to the Main script and adjust the parameters for the
# randomForest function. For detailed information on setup and additional
# options, refer to the Random Forest documentation:
# https://cran.r-project.org/web/packages/randomForest/index.html

################################################################################
#####                        User-defined parameters                       #####
################################################################################

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX274800"

# Optionally, use differential expression results to get genes that pass a
# minimum expression threshold. To use differential expression results, define
# the plot ID for the differential expression analysis you want to use.
# Filtering is recommended to reduce bias towards lowly abundant genes.
plot_id <- "52eaa19d-96bb-4b15-8793-a8283a93b759" # NULL skips this filtering

# Define the sample data column name to use for random forest.
sample_var <- "ajcc_pathologic_tumor_stage"

# Define the sample categories to use for random forest.
# These should be values from the sample_var column in your sample data.
sample_cat <- c("Stage IA", "Stage IIIA")

# Define a file path for the analysis plot (decision tree plot).
display_file_path <- paste0(experiment_id, "_decision_tree.png")

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(experiment_id, "_ranked_biomarkers.csv")

# Define a method to sort features by importance.
# Options: MeanDecreaseAccuracy or MeanDecreaseGini
sort_method <- "MeanDecreaseAccuracy"

# Define the number of important features to filter prior to the Classification
# and Regression Tree (CART).
n_features <- 50

# Include methods to describe your analysis.
plot_methods <- paste0(
    "Random Forest analysis was performed to identify important biomarkers ",
    "using the `randomForest` R package."
)

################################################################################
#####                            Main script                              #####
################################################################################

# Load required libraries
library(pluto)
library(edgeR)
library(randomForest)
library(dplyr)
library(rpart)
library(rpart.plot)

# Log into Pluto
pluto_login(api_token)

# Read in experiment assay data
assay_data <- pluto_read_assay_data(experiment_id)
rownames(assay_data) <- assay_data$gene_id
assay_data <- assay_data[, -1]

# Log2 CPM normalization of RNA-seq counts
assay_data_cpm <- as.data.frame(cpm(assay_data, log = FALSE))
assay_data_cpm <- round(log2(assay_data_cpm + 1), 3)

# Optionally filter assay_data by genes tested for differential expression
# Note: this does not use only differentially expressed genes
if (!is.null(plot_id)) {
    tryCatch(
        {
            dge_results <- pluto_read_results(
                experiment_id = experiment_id, plot_id = plot_id
            )
        },
        error = function(e) {
            stop(
                "Failed to read results for experiment ID: ", experiment_id,
                " and plot ID: ", plot_id
            )
        }
    )
    assay_data_cpm_filtered <- assay_data_cpm[
        rownames(assay_data_cpm) %in% dge_results$Gene_Symbol,
    ]
}

# Read in experiment sample data
sample_data <- pluto_read_sample_data(experiment_id)

# Filter sample data based on provided sample variable and categories
variable_data <- sample_data %>%
    filter(!!sym(sample_var) %in% sample_cat)

# Get the sample ids for the filtered data
sample_ids <- variable_data$sample_id

# Subset assay_data to only include columns corresponding to filtered sample_ids
assay_data_filtered <- assay_data_cpm_filtered[
    , colnames(assay_data_cpm_filtered) %in% sample_ids
]

# Transpose assay_data_filtered so that samples are rows and genes are columns
assay_data_filtered <- t(assay_data_filtered)

# Prepare variable for Random Forest (it should be a factor)
variable_vec <- factor(variable_data[[sample_var]])

# Run Random Forest for biomarker selection
set.seed(42) # For reproducibility
rf_model <- randomForest(
    assay_data_filtered,
    y = variable_vec,
    importance = TRUE,
    ntree = 500
)

# View the model summary
# print(rf_model)
# plot(rf_model)

# Feature importance (biomarker selection)
# Get the feature importance
importance_scores <- rf_model$importance

# Sort features by importance
sorted_importance <- importance_scores[
    order(importance_scores[, sort_method], decreasing = TRUE),
]

# Print top important biomarkers
print("Top biomarkers based on importance:")
head(sorted_importance, 10) # Top 10 important features

ranked_biomarkers <- as.data.frame(sorted_importance)
ranked_biomarkers <- cbind(
    gene_id = rownames(ranked_biomarkers), ranked_biomarkers
)
write.csv(ranked_biomarkers, results_file_path, row.names = FALSE, na = "")

# A single classification and regression tree (CART)
# Subset data to n_features features before classification
assay_data_cart <- data.frame(
    assay_data_filtered[, rownames(sorted_importance)[1:n_features]]
)
assay_data_cart$Group <- variable_vec # Directly assign the factor as Group
cart_model <- rpart(Group ~ .,
    data = assay_data_cart,
    method = "class"
)
# Plot the tree
png(display_file_path, width = 800, height = 600)
rpart.plot(cart_model, type = 5, extra = 101)
dev.off()

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)
