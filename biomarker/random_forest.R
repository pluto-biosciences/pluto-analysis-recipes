################################################################################
#####                 Random forest for biomarker selection                #####
################################################################################

# This script runs a Random Forest model to perform biomarker selection based
# on assay data and sample data from a Pluto Bio experiment.

# To use this script, simply update the user-defined parameters to reflect your
# specific experiment and analyses. You should not need to change anything else.

# For more advanced plot customization, scroll down to the Main script and edit
# parameters for the ggplot function. For more information, see the usage docs:
# https://ggplot2.tidyverse.org/

################################################################################
#####                        User-defined parameters                       #####
################################################################################

# !! Replace the following parameters with values relevant to your experiment !!

# You will need an API token to programmatically access your data in Pluto.
# For info, visit: https://help.pluto.bio/en/articles/creating-your-api-token
api_token <- "YOUR_API_TOKEN"

# The experiment ID can be found in the Workflow tab of your experiment in Pluto
# or in the URL, and always starts with "PLX".
experiment_id <- "PLX074620"

# Define the sample data column name to use for random forest.
sample_var <- "ajcc_pathologic_tumor_stage"

# Define the sample categories to use for random forest.
# These should be values from the sample_var column in your sample data.
sample_cat <- c("Stage IA", "Stage IB")

# Define a file path for the analysis plot (decision tree plot).
display_file_path <- paste0(experiment_id, "_decision_tree.png")

# Define a file path for the analysis results (CSV file).
results_file_path <- paste0(experiment_id, "_ranked_biomarkers.csv")

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
library(randomForest)
library(ggplot2)
library(dplyr)
library(ggraph)
library(igraph)

# Log into Pluto
pluto_login(api_token)

# Read in experiment assay data
assay_data <- pluto_read_assay_data(experiment_id)
rownames(assay_data) <- assay_data$gene_id
assay_data <- assay_data[, -1]

# Read in experiment sample data
sample_data <- pluto_read_sample_data(experiment_id)

# Filter sample data based on provided sample variable and categories
variable_data <- sample_data %>%
    filter(!!sym(sample_var) %in% sample_cat)

# Get the sample ids for the filtered data
sample_ids <- variable_data$sample_id

# Subset assay_data to only include columns corresponding to filtered sample_ids
assay_data_filtered <- assay_data[, colnames(assay_data) %in% sample_ids]

# Transpose assay_data_filtered so that samples are rows and genes are columns
assay_data_filtered <- t(assay_data_filtered)

# Prepare variable for Random Forest (it should be a factor)
variable_vec <- factor(variable_data[[sample_var]])

# Run Random Forest for biomarker selection
set.seed(42) # For reproducibility
rf_model <- randomForest(
    assay_data_filtered,
    y = variable_vec,
    importance = TRUE
)

# View the model summary
# print(rf_model)

# Plot the decision tree (take the first tree for simplicity)
# You can extract and plot any tree in the forest
tree <- getTree(rf_model, k = 1, labelVar = TRUE)
tree <- tree %>%
    tibble::rownames_to_column() %>%
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))
graph_frame <- data.frame(
    from = rep(tree$rowname, 2),
    to = c(tree$`left daughter`, tree$`right daughter`)
)
graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")

# Set node labels
V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
V(graph)$leaf_label <- as.character(tree$prediction)
V(graph)$split <- as.character(round(tree$`split point`, digits = 2))

# Plot
decision_tree_plot <- ggraph(graph, "dendrogram") +
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(
        aes(label = node_label),
        na.rm = TRUE, repel = TRUE
    ) +
    geom_node_label(
        aes(label = split),
        vjust = 2.5, na.rm = TRUE, fill = "white"
    ) +
    geom_node_label(
        aes(label = leaf_label, fill = leaf_label),
        na.rm = TRUE,
        repel = TRUE, colour = "white", fontface = "bold", show.legend = FALSE
    ) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 18)
    )

# Save decision tree plot
ggsave(
    display_file_path,
    height = 6, width = 12
)

# Feature importance (biomarker selection)
# Get the feature importance
importance_scores <- rf_model$importance

# Sort features by importance
sorted_importance <- importance_scores[
    order(importance_scores[, 1], decreasing = TRUE),
]

# Print top important biomarkers
print("Top biomarkers based on importance:")
head(sorted_importance, 10) # Top 10 important features

ranked_biomarkers <- as.data.frame(sorted_importance)
write.csv(ranked_biomarkers, results_file_path, row.names = FALSE, na = "")

# Push results back to Pluto
pluto_add_experiment_plot(
    experiment_id = experiment_id,
    display_file_path = display_file_path,
    analysis_name = analysis_name,
    plot_methods = plot_methods
)
