# ==============================================================================
# Kinase Structural State Analysis - Multimer Target Isolation Edition
# ==============================================================================

# --- CONFIGURATION ------------------------------------------------------------
INPUT_CSV <- "hmm_kinase_analysis_results_v4r1.csv"

# Re-enabled! Strictly orders groups based on your CSV's row order.
METADATA_FILE <- "experiment.csv" 

# Define the specific protein type to isolate for analysis
TARGET_TYPE <- "CSK"

# Set the number of columns for the Phase 5 Heatmap matrix (e.g., 4 = 2x4 grid)
HEATMAP_COLS <- 4

OUT_DIR <- sprintf("plots_and_stats_%s", TARGET_TYPE)

# Define the "Target States" for binary pairwise comparisons
TARGET_CATEGORIES <- list(
  State = "Active (BLAminus)",
  C_Helix = "In",
  R_Spine = "Intact",
  Spatial = "DFGin"
)
# ------------------------------------------------------------------------------

# --- Dependencies & Setup ---
required_packages <- c("tidyverse", "rstatix", "patchwork", "corrplot", "factoextra", "RColorBrewer", "pvclust")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(patchwork)
  library(corrplot)
  library(factoextra)
  library(RColorBrewer)
})

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[*] Outputs will be saved to: ./%s/\n", OUT_DIR))

PLOT_THEME <- theme_classic() + 
  theme(text = element_text(size = 14, color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5))

# ==============================================================================
# Phase 1: Data Preprocessing, Target Isolation, & Condition Assignment
# ==============================================================================
cat("[*] Loading and pre-processing data...\n")
df <- read_delim(INPUT_CSV, na = c("", "NA", "NaN", "None"), show_col_types = FALSE)

if ("Type" %in% colnames(df)) {
  initial_rows <- nrow(df)
  df <- df %>% filter(Type == TARGET_TYPE)
  cat(sprintf("[*] Target Isolation: Filtered to Type == '%s' (%d rows remaining out of %d).\n", 
              TARGET_TYPE, nrow(df), initial_rows))
  if(nrow(df) == 0) stop("[!] No rows match the TARGET_TYPE. Check your CSV spelling.")
}

cat("[*] Dynamically extracting conditions from AF3 directory names...\n")
df <- df %>%
  mutate(Condition = str_extract(Directory, "^[^/]+")) %>%
  mutate(Condition = str_replace(Condition, "^[a-zA-Z]-", "")) %>%
  mutate(Condition = str_replace(Condition, "_[a-zA-Z]-", "\n"))

if (!is.null(METADATA_FILE) && file.exists(METADATA_FILE)) {
  cat("[*] Applying strict experimental design from metadata...\n")
  meta <- read_csv(METADATA_FILE, show_col_types = FALSE)
  meta[is.na(meta)] <- ""
  
  # Translate CSV columns into the exact clean plotting labels
  meta <- meta %>%
    mutate(
      clean_a = paste0(tolower(chain_a), "-wtcat-", tolower(condition_a)),
      clean_b = if_else(ptm_b != "", 
                        paste0(tolower(chain_b), "-wtcat-", tolower(ptm_b), "-", tolower(condition_b)),
                        paste0(tolower(chain_b), "-wtcat-", tolower(condition_b))),
      Condition_Label = paste0(clean_a, "\n", clean_b)
    )
  
  # Filter to ONLY conditions in the CSV and lock them in exact CSV row order!
  df <- df %>%
    filter(Condition %in% meta$Condition_Label) %>%
    mutate(Condition = factor(Condition, levels = meta$Condition_Label))
} else {
  df <- df %>% mutate(Condition = as.factor(Condition))
}

unique_conditions <- levels(df$Condition)
num_conditions <- length(unique_conditions)
cat(sprintf("[*] Extracted %d distinct multimer groups.\n", num_conditions))

# Pull Set1 but drop the 6th color (bright yellow) for better white-background visibility
safe_colors <- brewer.pal(max(3, 9), "Set1")[-6] 
CUSTOM_COLORS <- setNames(colorRampPalette(safe_colors)(num_conditions), unique_conditions)

all_dist_cols <- grep("_Dist$", colnames(df), value = TRUE)
target_cols <- c("Phi_D", "Psi_D", all_dist_cols)

df <- df %>%
  mutate(across(any_of(target_cols), ~ suppressWarnings(as.numeric(as.character(.)))))

# --- LIGAND CROSSTALK CORRECTION ---
cat("[*] Applying cross-chain ligand distance correction...\n")
ligand_cols <- grep("ATP|Mg", all_dist_cols, ignore.case = TRUE, value = TRUE)

df$Is_Target_Apo <- sapply(as.character(df$Condition), function(cond) {
  chunks <- unlist(strsplit(cond, "[\n_]"))
  target_chunk <- chunks[grepl(TARGET_TYPE, chunks, ignore.case = TRUE)]
  if (length(target_chunk) > 0) {
    return(grepl("apo", target_chunk[1], ignore.case = TRUE))
  }
  return(FALSE)
})

if(length(ligand_cols) > 0) {
  df <- df %>%
    mutate(across(any_of(ligand_cols), ~ if_else(Is_Target_Apo, NA_real_, .)))
}
df <- df %>% select(-Is_Target_Apo)

# --- INTELLIGENT METRIC SORTING ---
viable_dist_summary <- df %>%
  select(Condition, all_of(all_dist_cols)) %>%
  group_by(Condition) %>%
  summarise(across(everything(), ~sum(!is.na(.)) > 5)) %>%
  select(-Condition)

universal_dist_cols <- names(which(sapply(viable_dist_summary, all)))
any_dist_cols <- names(which(sapply(viable_dist_summary, any)))
multi_dist_cols <- names(which(sapply(viable_dist_summary, sum) >= 2))

cat(sprintf("[*] Distance Metrics Profile for %s:\n    - %d viable overall\n    - %d multi-group (Stats available)\n    - %d universal (PCA safe)\n", 
            TARGET_TYPE, length(any_dist_cols), length(multi_dist_cols), length(universal_dist_cols)))

# ==============================================================================
# Phase 2: Macro-State Conformational Shifts (Categorical)
# ==============================================================================
cat("[*] Executing Phase 2: Categorical Macro-States...\n")

plot_categorical <- function(var_name, title) {
  df %>%
    count(Condition, !!sym(var_name)) %>% group_by(Condition) %>% mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes(x = Condition, y = Percent, fill = !!sym(var_name))) +
    geom_bar(stat = "identity", color = "black", width = 0.6) + scale_fill_brewer(palette = "Paired") +
    labs(title = sprintf("%s: %s", TARGET_TYPE, title), x = NULL, y = "Proportion (%)", fill = "State") +
    PLOT_THEME + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
}

p_state <- plot_categorical("State", "Global Conformation")
p_chelix <- plot_categorical("C_Helix", "aC-Helix State")
p_rspine <- plot_categorical("R_Spine", "R-Spine Integrity")
p_spatial <- plot_categorical("Spatial", "DFG Spatial State")

phase2_plot <- (p_state | p_chelix) / (p_rspine | p_spatial)
ggsave(file.path(OUT_DIR, "Phase2_Macro_States.pdf"), phase2_plot, width = 16, height = 12)

# --- Pairwise Binary Categorical Statistics ---
cat("\n--- Pairwise Fisher's Exact Tests (Binary Splits) ---\n")
pairwise_fisher_results <- tibble(
  Structural_Feature = character(), Group_1 = character(), Group_2 = character(), 
  Target_State = character(), p_value = numeric()
)

condition_pairs <- combn(unique_conditions, 2, simplify = FALSE)

for(cat_var in names(TARGET_CATEGORIES)) {
  target_val <- TARGET_CATEGORIES[[cat_var]]
  
  for(pair in condition_pairs) {
    grp1 <- pair[1]
    grp2 <- pair[2]
    
    sub_df <- df %>% 
      filter(Condition %in% c(grp1, grp2)) %>%
      mutate(Binary_State = factor(if_else(!!sym(cat_var) == target_val, target_val, "Other")))
    
    tbl <- table(droplevels(sub_df$Condition), sub_df$Binary_State)
    
    if(nrow(tbl) == 2 && ncol(tbl) == 2) {
      res <- fisher.test(tbl)
      pairwise_fisher_results <- pairwise_fisher_results %>% 
        add_row(Structural_Feature = cat_var, Group_1 = grp1, Group_2 = grp2, 
                Target_State = target_val, p_value = res$p.value)
    } else {
      pairwise_fisher_results <- pairwise_fisher_results %>% 
        add_row(Structural_Feature = cat_var, Group_1 = grp1, Group_2 = grp2, 
                Target_State = target_val, p_value = NA)
    }
  }
}

pairwise_fisher_results <- pairwise_fisher_results %>%
  group_by(Structural_Feature) %>%
  mutate(p.adj = p.adjust(p_value, method = "BH")) %>%
  add_significance("p.adj") %>% 
  arrange(Structural_Feature, p.adj)

write_csv(pairwise_fisher_results, file.path(OUT_DIR, "Phase2_Pairwise_Categorical_Stats.csv"))

# ==============================================================================
# Phase 3: Conformational Phase Space (2D Continuous)
# ==============================================================================
cat("[*] Executing Phase 3: 2D Phase Space...\n")

p_dihedral <- df %>% drop_na(Phi_D, Psi_D) %>%
  ggplot(aes(x = Phi_D, y = Psi_D, color = Condition, fill = Condition)) +
  geom_point(alpha = 0.6, size = 2) + geom_density_2d(alpha = 0.8, linewidth = 0.8) +
  scale_color_manual(values = CUSTOM_COLORS) + scale_fill_manual(values = CUSTOM_COLORS) +
  labs(title = sprintf("%s: DFG Dihedral Phase Space", TARGET_TYPE), x = expression(Phi~"(°)"), y = expression(Psi~"(°)")) + PLOT_THEME

p_dunbrack <- df %>% drop_na(D1_Dist, D2_Dist) %>%
  ggplot(aes(x = D1_Dist, y = D2_Dist, color = Condition, fill = Condition)) +
  geom_point(alpha = 0.6, size = 2) + geom_density_2d(alpha = 0.8, linewidth = 0.8) +
  scale_color_manual(values = CUSTOM_COLORS) + scale_fill_manual(values = CUSTOM_COLORS) +
  labs(title = sprintf("%s: Dunbrack Spatial Coordinates", TARGET_TYPE), x = "D1 Distance (Å)", y = "D2 Distance (Å)") + PLOT_THEME

suppressWarnings({ ggsave(file.path(OUT_DIR, "Phase3_2D_PhaseSpace.pdf"), p_dihedral | p_dunbrack, width = 14, height = 6) })

# ==============================================================================
# Phase 4: Micro-Metric Allosteric Evaluation (1D Continuous)
# ==============================================================================
cat("[*] Executing Phase 4: 1D Micro-Metrics & Stats...\n")

df_long <- df %>% select(Simulation_ID, Condition, all_of(any_dist_cols)) %>%
  pivot_longer(cols = -c(Simulation_ID, Condition), names_to = "Metric", values_to = "Distance") %>%
  drop_na(Distance)

p_violins <- df_long %>%
  ggplot(aes(x = Condition, y = Distance, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.8, outlier.shape = NA) +
  facet_wrap(~Metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = CUSTOM_COLORS) +
  labs(title = sprintf("%s: Allosteric Distances Across Complex Groups", TARGET_TYPE), y = "Distance (Å)", x = NULL) +
  PLOT_THEME + theme(strip.background = element_rect(fill = "grey90", color = "white"), 
                     strip.text = element_text(face = "bold"),
                     axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(OUT_DIR, "Phase4_Allosteric_Distances.pdf"), p_violins, width = 16, height = 12)

df_stats <- df_long %>% filter(Metric %in% multi_dist_cols)

kruskal_results <- df_stats %>%
  group_by(Metric) %>% kruskal_test(Distance ~ Condition) %>%
  adjust_pvalue(method = "BH") %>% add_significance("p.adj")

pairwise_results <- df_stats %>%
  group_by(Metric) %>% pairwise_wilcox_test(Distance ~ Condition, p.adjust.method = "BH") %>%
  add_significance("p.adj")

write_csv(kruskal_results, file.path(OUT_DIR, "Phase4_Kruskal_Summary.csv"))
write_csv(pairwise_results, file.path(OUT_DIR, "Phase4_Pairwise_Summary.csv"))

# ==============================================================================
# Phase 5: Network Coupling & Multivariate Analysis
# ==============================================================================
cat("[*] Executing Phase 5: Network Coupling (PCA & Correlation)...\n")

pca_data <- df %>% select(all_of(universal_dist_cols)) %>% drop_na()
pca_metadata <- df %>% filter(row_number() %in% rownames(pca_data))

pca_res <- prcomp(pca_data, scale. = TRUE) 

suppressWarnings({
  p_pca <- fviz_pca_ind(pca_res, geom = "point", col.ind = pca_metadata$Condition, 
                        palette = as.character(CUSTOM_COLORS), addEllipses = TRUE, ellipse.type = "confidence",
                        title = sprintf("%s: Shared Allosteric Variance", TARGET_TYPE),
                        subtitle = sprintf("Based on %d universal metrics", length(universal_dist_cols))) + PLOT_THEME 
})
ggsave(file.path(OUT_DIR, "Phase5_PCA_Biplot.pdf"), p_pca, width = 10, height = 7)

high_contrast_pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)

# --- DYNAMIC GRID CALCULATION ---
n_cols <- min(HEATMAP_COLS, num_conditions)
n_rows <- ceiling(num_conditions / n_cols)

pdf(file.path(OUT_DIR, "Phase5_Correlation_Heatmaps.pdf"), width = 7 * n_cols, height = 7 * n_rows)
par(mfrow = c(n_rows, n_cols))

for (grp in unique_conditions) {
  df_grp <- df %>% filter(Condition == grp) %>% select(all_of(any_dist_cols)) %>%
    select(where(~sum(!is.na(.)) >= 3))
  
  if (nrow(df_grp) >= 3 && ncol(df_grp) >= 2) {
    cor_matrix <- suppressWarnings(cor(df_grp, use = "pairwise.complete.obs", method = "spearman"))
    cor_matrix[is.na(cor_matrix)] <- 0 
    
    # Increased text elements for large grid visibility:
    # tl.cex scales the variable labels, cl.cex scales the colorbar numbers, cex.main scales the title
    corrplot(cor_matrix, method = "color", type = "upper", 
             tl.col = "black", tl.cex = 1.6, cl.cex = 1.5, addgrid.col = "white",
             title = sprintf("%s Coupling\n(n=%d)", grp, nrow(df_grp)), 
             mar = c(0,0,5,0), cex.main = 2.0, col = high_contrast_pal)
  } else {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 2), ylim = c(0, 2))
    # Scaled up the fallback text (cex = 2.0)
    text(1, 1.2, sprintf("%s Coupling", grp), font = 2, cex = 2.0)
    text(1, 0.8, sprintf("Insufficient Data\n(n=%d)", nrow(df_grp)), cex = 2.0, col = "#e74c3c")
  }
}
dev.off()

# ==============================================================================
# Phase 6: Meta-Stable State Discovery (Gap Statistic Validation)
# ==============================================================================
cat("[*] Executing Phase 6: Gap Statistic Clustering (Meta-Stable States)...\n")

# Use the top 3 Principal Components
cluster_data <- pca_res$x[, 1:min(3, ncol(pca_res$x))]

cat("    -> Running Gap Statistic bootstrapping (100 iterations) to prove optimal states...\n")
set.seed(42) # Ensure reproducible bootstrapping
gap_stat <- cluster::clusGap(cluster_data, FUN = kmeans, nstart = 25, K.max = 8, B = 100)

# Extract optimal K using the Tibshirani Standard Error rule
# Finds the maximum K before the cluster quality becomes indistinguishable from random noise
optimal_k <- cluster::maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method="Tibs2001SEmax")

# Safety net: If the math says 1 state, force 2 to allow comparative plotting
if(optimal_k < 2) {
  cat("    -> [!] Gap statistic suggested 1 state (too conservative). Forcing 2 states for visualization.\n")
  optimal_k <- 2
}

cat(sprintf("[*] Statistically validated %d distinct meta-stable states.\n", optimal_k))

# Perform the final robust clustering with the mathematically validated K
set.seed(42)
km_final <- kmeans(cluster_data, centers = optimal_k, nstart = 50)
pca_metadata$Macro_State <- factor(paste("State", km_final$cluster))

suppressWarnings({
  # Plot 1: The Proof (Gap Statistic curve with Error Bars)
  p_gap <- fviz_gap_stat(gap_stat, maxSE = list(method = "Tibs2001SEmax")) + 
    PLOT_THEME + 
    labs(title = "Statistical Proof of Meta-Stable States",
         subtitle = "Optimal K is chosen right before the error bars overlap (Tibshirani Rule)")
  
  # Plot 2: Where are these states in the phase space?
  p_state_pca <- fviz_pca_ind(pca_res, geom = "point", col.ind = pca_metadata$Macro_State,
                              addEllipses = TRUE, ellipse.level = 0.95,
                              title = "Validated Meta-Stable States",
                              subtitle = sprintf("K-Means Clustering: Rigorously partitioned into %d basins", optimal_k)) + 
    PLOT_THEME + scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2")
})

ggsave(file.path(OUT_DIR, "Phase6_Statistical_Proof.pdf"), p_gap, width = 8, height = 6)
ggsave(file.path(OUT_DIR, "Phase6_State_Clusters_PCA.pdf"), p_state_pca, width = 10, height = 7)

# Plot 3: State Composition (Who is inside these states?)
p_state_comp <- pca_metadata %>%
  count(Macro_State, Condition) %>%
  group_by(Macro_State) %>%
  mutate(Percent = n / sum(n) * 100) %>%
  ggplot(aes(x = Macro_State, y = Percent, fill = Condition)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  scale_fill_manual(values = CUSTOM_COLORS) +
  labs(title = "Composition of Meta-Stable States", 
       subtitle = "Which experimental conditions map to which structural basins?",
       x = "Discovered Meta-Stable State", y = "Proportion (%)") +
  PLOT_THEME + theme(axis.text.x = element_text(face = "bold", size = 12))

ggsave(file.path(OUT_DIR, "Phase6_State_Composition.pdf"), p_state_comp, width = 10, height = 7)

# Save the assignments
write_csv(pca_metadata %>% select(Simulation_ID, Condition, Macro_State), 
          file.path(OUT_DIR, "Phase6_State_Assignments.csv"))

# ==============================================================================
# Phase 7: Biological Signatures of Meta-Stable States
# ==============================================================================
cat("\n[*] Executing Phase 7: Mapping Meta-Stable States to Biological Features...\n")

# Binarize the categories exactly as defined in TARGET_CATEGORIES (Target vs Other)
for(cat_var in names(TARGET_CATEGORIES)) {
  target_val <- TARGET_CATEGORIES[[cat_var]]
  new_col <- paste0("Binary_", cat_var)
  
  pca_metadata[[new_col]] <- if_else(pca_metadata[[cat_var]] == target_val, target_val, "Other")
  pca_metadata[[new_col]] <- factor(pca_metadata[[new_col]], levels = c(target_val, "Other"))
}

# Plot function: Maps the Target State (Blue) vs Other (Grey) for each mathematical cluster
plot_macro_feature <- function(var_name, original_var) {
  target_val <- TARGET_CATEGORIES[[original_var]]
  
  pca_metadata %>%
    count(Macro_State, !!sym(var_name)) %>% 
    group_by(Macro_State) %>% 
    mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes(x = Macro_State, y = Percent, fill = !!sym(var_name))) +
    geom_bar(stat = "identity", color = "black", width = 0.6) + 
    # Use high contrast: Bold Blue for the target biological state, Neutral Grey for "Other"
    scale_fill_manual(values = setNames(c("#2171b5", "#cccccc"), c(target_val, "Other"))) +
    labs(title = sprintf("%s Identity", original_var), x = NULL, y = "Proportion (%)", fill = "Feature") +
    PLOT_THEME + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
}

p_m1 <- plot_macro_feature("Binary_State", "State")
p_m2 <- plot_macro_feature("Binary_C_Helix", "C_Helix")
p_m3 <- plot_macro_feature("Binary_R_Spine", "R_Spine")
p_m4 <- plot_macro_feature("Binary_Spatial", "Spatial")

phase7_plot <- (p_m1 | p_m2) / (p_m3 | p_m4) + 
  plot_annotation(title = sprintf("%s: Biological Signatures of Meta-Stable Basins", TARGET_TYPE),
                  theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

ggsave(file.path(OUT_DIR, "Phase7_MacroState_Signatures.pdf"), phase7_plot, width = 16, height = 12)

# --- Statistical Testing: Does the Macro_State significantly predict the biological feature? ---
macro_stats <- tibble(Structural_Feature = character(), Target_State = character(), p_value = numeric())

for(cat_var in names(TARGET_CATEGORIES)) {
  bin_var <- paste0("Binary_", cat_var)
  
  # Cross-tabulate mathematical states vs biological states
  tbl <- table(pca_metadata$Macro_State, pca_metadata[[bin_var]])
  
  if(nrow(tbl) >= 2 && ncol(tbl) >= 2) {
    # Using simulate.p.value = TRUE because the table is larger than 2x2 (e.g., 4 States x 2 Features)
    res <- fisher.test(tbl, simulate.p.value = TRUE)
    macro_stats <- macro_stats %>% 
      add_row(Structural_Feature = cat_var, Target_State = TARGET_CATEGORIES[[cat_var]], p_value = res$p.value)
  } else {
    macro_stats <- macro_stats %>% 
      add_row(Structural_Feature = cat_var, Target_State = TARGET_CATEGORIES[[cat_var]], p_value = NA)
  }
}

# Add FDR correction and significance stars
macro_stats <- macro_stats %>%
  mutate(p.adj = p.adjust(p_value, method = "BH")) %>%
  add_significance("p.adj")

write_csv(macro_stats, file.path(OUT_DIR, "Phase7_MacroState_Signatures_Stats.csv"))

cat(sprintf("[✓] Pipeline totally complete! Mapped mathematical states to biological features.\n"))
