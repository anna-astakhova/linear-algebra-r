library(ggplot2)

# --- Plot 1: scree plot ---

pca_result <- compute_pca(as.matrix(iris[, 1:4]))

# scree plot
scree_df <- data.frame(
  PC       = 1:4,
  Variance = pca_result$var_explained * 100,
  Cumulative = pca_result$cum_var_explained * 100
)

ggplot(scree_df, aes(x = PC)) +
  geom_col(aes(y = Variance), fill = "#7F77DD", alpha = 0.8, width = 0.6) +
  geom_line(aes(y = Cumulative), colour = "#1D9E75", linewidth = 1) +
  geom_point(aes(y = Cumulative), colour = "#1D9E75", size = 3) +
  geom_hline(yintercept = 95, linetype = "dashed",
             colour = "#BA7517", alpha = 0.7) +
  annotate("text", x = 3.5, y = 96.5, label = "95% threshold",
           size = 3, colour = "#BA7517") +
  scale_y_continuous(limits = c(0, 105)) +
  labs(
    title    = "Scree Plot: Variance Explained by Each PC",
    subtitle = "Bars = individual; Line = cumulative",
    x        = "Principal Component",
    y        = "Variance Explained (%)"
  ) +
  theme_minimal(base_size = 12)
ggsave("figures/pca_scree_plot1.png", dpi=150)

# --- Plot 2: score plot coloured by species ---

# score plot coloured by species
scores_df <- data.frame(
  PC1     = pca_result$scores[, 1],
  PC2     = pca_result$scores[, 2],
  Species = iris$Species
)

# compute % variance for axis labels
pc1_var <- round(pca_result$var_explained[1] * 100, 1)
pc2_var <- round(pca_result$var_explained[2] * 100, 1)

ggplot(scores_df, aes(x = PC1, y = PC2, colour = Species)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = c("#7F77DD", "#1D9E75", "#BA7517")) +
  labs(
    title    = "PCA of Iris Dataset — PC1 vs PC2",
    subtitle = "Implemented from scratch using eigendecomposition of the covariance matrix",
    x        = paste0("PC1 (", pc1_var, "% variance)"),
    y        = paste0("PC2 (", pc2_var, "% variance)")
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

ggsave("figures/pca_iris_biplot.png", width=7, height=5, dpi=150)

# --- Plot 3: loadings bar chart ---
var_names <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")

loadings_df <- data.frame(
  Variable  = rep(var_names, 2),
  Loading   = c(pca_result$loadings[, 1], pca_result$loadings[, 2]),
  Component = rep(c("PC1", "PC2"), each = 4)
)

ggplot(loadings_df, aes(x = Variable, y = Loading, fill = Component)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.6) +
  scale_fill_manual(values = c("#7F77DD", "#1D9E75")) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.3) +
  labs(
    title    = "PCA Loadings: Contribution of Each Variable",
    subtitle = "Large absolute loading = strong contribution to that PC",
    x        = NULL,
    y        = "Loading"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave("figures/pca_loadings_bar_chart.png", dpi=150)


#-------------------------------------------------------
# --- Plot 1: scree plot ---
scree_df <- data.frame(
  PC         = 1:10,
  Variance   = pca_result$var_explained[1:10] * 100,
  Cumulative = pca_result$cum_var_explained[1:10] * 100
)

ggplot(scree_df, aes(x = PC)) +
  geom_col(aes(y = Variance), fill = "#7F77DD", alpha = 0.8, width = 0.6) +
  geom_line(aes(y = Cumulative), colour = "#1D9E75", linewidth = 1) +
  geom_point(aes(y = Cumulative), colour = "#1D9E75", size = 3) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "#BA7517") +
  labs(
    title    = "Neural Population PCA — Variance Explained",
    subtitle = sprintf("%d neurons, %d trials", ncol(neural_mat), nrow(neural_mat)),
    x        = "Principal Component",
    y        = "Variance Explained (%)"
  ) +
  theme_minimal(base_size = 12)

ggsave("figures/neural_pca_scree.png", width = 7, height = 4, dpi = 150)

# --- Plot 2: PC1 vs PC2 coloured by condition ---

scores_df <- data.frame(
  PC1 = pca_result$scores[, 1],
  PC2 = pca_result$scores[, 2],
  PC3 = pca_result$scores[, 3]
)

pc1_var <- round(pca_result$var_explained[1] * 100, 1)
pc2_var <- round(pca_result$var_explained[2] * 100, 1)

ggplot(scores_df, aes(x = PC1, y = PC2, colour = PC3)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_gradient2(low = "#7F77DD", mid = "grey80",
                         high = "#1D9E75", midpoint = 0) +
  labs(
    title    = "Neural Population Trajectory — PC1 vs PC2",
    subtitle = "Each point = one trial; colour = PC3 score",
    x        = paste0("PC1 (", pc1_var, "% variance)"),
    y        = paste0("PC2 (", pc2_var, "% variance)"),
    colour   = "PC3"
  ) +
  theme_minimal(base_size = 12)

ggsave("figures/neural_pca_scores.png", width = 7, height = 5, dpi = 150)

# --- Plot 3: loadings — which neurons drive each PC? ---
# This is the biological interpretation plot

n_neurons_show <- min(20, ncol(neural_mat))  # show top 20

loadings_df <- data.frame(
  Neuron  = rep(paste0("N", 1:n_neurons_show), 2),
  Loading = c(pca_result$loadings[1:n_neurons_show, 1],
              pca_result$loadings[1:n_neurons_show, 2]),
  PC      = rep(c("PC1", "PC2"), each = n_neurons_show)
)

ggplot(loadings_df, aes(x = Neuron, y = Loading, fill = PC)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c("#7F77DD", "#1D9E75")) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  labs(
    title    = "PCA Loadings: Neuron Contributions to PC1 and PC2",
    subtitle = "Large absolute loading = neuron strongly drives that component",
    x        = "Neuron", y = "Loading"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/neural_pca_loadings.png", width = 8, height = 4, dpi = 150)
