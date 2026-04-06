# ==============================================================================
# Kinase Structural State Analysis
# ==============================================================================

# --- CONFIGURATION ------------------------------------------------------------
# 1. Main Input File
INPUT_CSV <- "hmm_kinase_analysis_results_v4r1.csv"

# 2. Metadata Input (Optional)
# Set to a file path (e.g., "metadata.tsv" or "design.csv") to merge external info.
# Set to NULL if you want to extract conditions directly from the INPUT_CSV strings.
METADATA_FILE <- NULL 

# 3. Contrast / Group Definitions
# If METADATA_FILE = NULL: The script searches this column in the INPUT_CSV for the Labels below.
# If METADATA_FILE = "file.csv": This must match the exact column name in the metadata file defining the groups.
CONDITION_COL <- "Directory" 

# Define the two comparison groups. 
# (If using regex extraction, these act as the case-insensitive search patterns).
GROUP_A_LABEL <- "Apo"
GROUP_B_LABEL <- "Holo"

OUT_DIR <- "plots_and_stats"
# ------------------------------------------------------------------------------

# --- Dependencies & Setup ---
required_packages <- c("tidyverse", "rstatix", "patchwork", "corrplot", "factoextra", "RColorBrewer")
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

# Create output directory if it doesn't exist
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[*] Outputs will be saved to: ./%s/\n", OUT_DIR))

CUSTOM_COLORS <- setNames(c("#3498db", "#e74c3c"), c(GROUP_A_LABEL, GROUP_B_LABEL))
PLOT_THEME <- theme_classic() + 
  theme(text = element_text(size = 14, color = "black"),
        axis.text = element_text(color = "black"),
        legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5))

# ==============================================================================
# Phase 1: Data Preprocessing & Condition Assignment
# ==============================================================================
cat("[*] Loading and pre-processing data...\n")
df <- read_csv(INPUT_CSV, na = c("", "NA", "NaN", "None"), show_col_types = FALSE)

if (!is.null(METADATA_FILE) && file.exists(METADATA_FILE)) {
  delim <- if(grepl("\\.tsv$", METADATA_FILE, ignore.case=TRUE)) "\t" else ","
  meta <- read_delim(METADATA_FILE, delim = delim, show_col_types = FALSE)
  join_key <- intersect(colnames(df), colnames(meta))[1]
  df <- df %>% left_join(meta, by = join_key)
  df <- df %>% rename(Condition = !!sym(CONDITION_COL)) %>%
    filter(Condition %in% c(GROUP_A_LABEL, GROUP_B_LABEL))
} else {
  cat(sprintf("[*] Extracting contrasts via regex on column: %s\n", CONDITION_COL))
  df <- df %>%
    mutate(Condition = case_when(
      str_detect(tolower(!!sym(CONDITION_COL)), tolower(GROUP_A_LABEL)) ~ GROUP_A_LABEL,
      str_detect(tolower(!!sym(CONDITION_COL)), tolower(GROUP_B_LABEL)) ~ GROUP_B_LABEL,
      TRUE ~ "Other"
    )) %>%
    filter(Condition %in% c(GROUP_A_LABEL, GROUP_B_LABEL))
}

df <- df %>% mutate(Condition = factor(Condition, levels = c(GROUP_A_LABEL, GROUP_B_LABEL)))

# Force numerical typing
all_dist_cols <- grep("_Dist$", colnames(df), value = TRUE)
target_cols <- c("Phi_D", "Psi_D", all_dist_cols)

df <- df %>%
  mutate(across(any_of(target_cols), ~ suppressWarnings(as.numeric(as.character(.)))))

viable_dist_cols <- df %>%
  select(Condition, all_of(all_dist_cols)) %>%
  group_by(Condition) %>%
  summarise(across(everything(), ~sum(!is.na(.)) > 5)) %>% 
  summarise(across(-Condition, all)) %>%
  unlist()

dist_cols <- names(viable_dist_cols[viable_dist_cols == TRUE])
cat(sprintf("[*] Auto-detected %d shared distance metrics for statistical comparison.\n", length(dist_cols)))

# ==============================================================================
# Phase 2: Macro-State Conformational Shifts (Categorical)
# ==============================================================================
cat("[*] Executing Phase 2: Categorical Macro-States...\n")

plot_categorical <- function(var_name, title) {
  df %>%
    count(Condition, !!sym(var_name)) %>%
    group_by(Condition) %>%
    mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes(x = Condition, y = Percent, fill = !!sym(var_name))) +
    geom_bar(stat = "identity", color = "black", width = 0.6) +
    scale_fill_brewer(palette = "Paired") +
    labs(title = title, x = NULL, y = "Proportion (%)", fill = "State") +
    PLOT_THEME + theme(legend.position = "right")
}

p_state <- plot_categorical("State", "Global Conformation")
p_chelix <- plot_categorical("C_Helix", "aC-Helix State")
p_rspine <- plot_categorical("R_Spine", "R-Spine Integrity")
p_spatial <- plot_categorical("Spatial", "DFG Spatial State")

phase2_plot <- (p_state | p_chelix) / (p_rspine | p_spatial)
ggsave(file.path(OUT_DIR, "Phase2_Macro_States.pdf"), phase2_plot, width = 12, height = 10)

cat("\n--- Fisher's Exact Tests (Categorical Shifts) ---\n")
fisher_results <- tibble(Structural_Feature = character(), p_value = numeric())

for(cat_var in c("State", "C_Helix", "R_Spine", "Spatial")) {
  tbl <- table(df$Condition, df[[cat_var]])
  res <- fisher.test(tbl, simulate.p.value = TRUE)
  
  fisher_results <- fisher_results %>% add_row(Structural_Feature = cat_var, p_value = res$p.value)
  cat(sprintf("%s vs Condition: p-value = %.4e\n", cat_var, res$p.value))
}

write_csv(fisher_results, file.path(OUT_DIR, "Phase2_Categorical_Stats.csv"))

# ==============================================================================
# Phase 3: Conformational Phase Space (2D Continuous)
# ==============================================================================
cat("\n[*] Executing Phase 3: 2D Phase Space...\n")

p_dihedral <- df %>% drop_na(Phi_D, Psi_D) %>%
  ggplot(aes(x = Phi_D, y = Psi_D, color = Condition, fill = Condition)) +
  geom_point(alpha = 0.6, size = 2) + geom_density_2d(alpha = 0.8, linewidth = 0.8) +
  scale_color_manual(values = CUSTOM_COLORS) + scale_fill_manual(values = CUSTOM_COLORS) +
  labs(title = "DFG Dihedral Phase Space", x = expression(Phi~"(°)"), y = expression(Psi~"(°)")) + PLOT_THEME

p_dunbrack <- df %>% drop_na(D1_Dist, D2_Dist) %>%
  ggplot(aes(x = D1_Dist, y = D2_Dist, color = Condition, fill = Condition)) +
  geom_point(alpha = 0.6, size = 2) + geom_density_2d(alpha = 0.8, linewidth = 0.8) +
  scale_color_manual(values = CUSTOM_COLORS) + scale_fill_manual(values = CUSTOM_COLORS) +
  labs(title = "Dunbrack Spatial Coordinates", x = "D1 Distance (Å)", y = "D2 Distance (Å)") + PLOT_THEME

suppressWarnings({
  ggsave(file.path(OUT_DIR, "Phase3_2D_PhaseSpace.pdf"), p_dihedral | p_dunbrack, width = 12, height = 5)
})

# ==============================================================================
# Phase 4: Micro-Metric Allosteric Evaluation (1D Continuous)
# ==============================================================================
cat("[*] Executing Phase 4: 1D Micro-Metrics & Stats...\n")

df_long <- df %>% select(Simulation_ID, Condition, all_of(dist_cols)) %>%
  pivot_longer(cols = -c(Simulation_ID, Condition), names_to = "Metric", values_to = "Distance") %>%
  mutate(Distance = suppressWarnings(as.numeric(as.character(Distance)))) %>%
  drop_na(Distance)

p_violins <- df_long %>%
  ggplot(aes(x = Condition, y = Distance, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.2, color = "black", alpha = 0.8, outlier.shape = NA) +
  facet_wrap(~Metric, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = CUSTOM_COLORS) +
  labs(title = sprintf("Allosteric Network Distances: %s vs %s", GROUP_A_LABEL, GROUP_B_LABEL), y = "Distance (Å)", x = NULL) +
  PLOT_THEME + theme(strip.background = element_rect(fill = "grey90", color = "white"), strip.text = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "Phase4_Allosteric_Distances.pdf"), p_violins, width = 14, height = 10)

stat_results <- df_long %>%
  group_by(Metric) %>%
  wilcox_test(Distance ~ Condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

final_stats <- stat_results

tryCatch({
  effect_sizes <- df_long %>%
    group_by(Metric) %>%
    wilcox_effsize(Distance ~ Condition) %>% 
    select(Metric, Effect_Size_r = effsize)
  
  final_stats <- final_stats %>% left_join(effect_sizes, by = "Metric")
}, error = function(e) {
  cat("[!] Note: Effect size calculation skipped due to zero-variance data limits. Saving p-values only.\n")
})

final_stats <- final_stats %>% arrange(p.adj)
write_csv(final_stats, file.path(OUT_DIR, "Phase4_Statistical_Summary.csv"))

# ==============================================================================
# Phase 5: Network Coupling & Multivariate Analysis
# ==============================================================================
cat("[*] Executing Phase 5: Network Coupling (PCA & Correlation)...\n")

pca_data <- df %>% select(all_of(dist_cols)) %>% drop_na()
pca_metadata <- df %>% filter(row_number() %in% rownames(pca_data))

pca_res <- prcomp(pca_data, scale. = TRUE) 
p_pca <- fviz_pca_ind(pca_res, geom = "point", col.ind = pca_metadata$Condition, 
                      palette = as.character(CUSTOM_COLORS), addEllipses = TRUE, ellipse.type = "confidence",
                      title = "PCA: Global Allosteric Variance") + 
  PLOT_THEME + 
  labs(caption = "Note: Circled areas (ellipses) represent the 95% confidence intervals of the group centroids.")

ggsave(file.path(OUT_DIR, "Phase5_PCA_Biplot.pdf"), p_pca, width = 7, height = 6)

df_a <- df %>% filter(Condition == GROUP_A_LABEL) %>% select(all_of(dist_cols)) %>% drop_na()
df_b <- df %>% filter(Condition == GROUP_B_LABEL) %>% select(all_of(dist_cols)) %>% drop_na()

# Create a high-contrast Divergent Palette (Deep Blue -> White -> Deep Red)
high_contrast_pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)

pdf(file.path(OUT_DIR, "Phase5_Correlation_Heatmaps.pdf"), width = 14, height = 7)
par(mfrow = c(1, 2))

corrplot(cor(df_a, method = "spearman"), method = "color", type = "upper", 
         tl.col = "black", tl.cex = 0.8, addgrid.col = "white",
         title = sprintf("%s Network Coupling", GROUP_A_LABEL), mar = c(0,0,2,0), 
         col = high_contrast_pal)

corrplot(cor(df_b, method = "spearman"), method = "color", type = "upper", 
         tl.col = "black", tl.cex = 0.8, addgrid.col = "white",
         title = sprintf("%s Network Coupling", GROUP_B_LABEL), mar = c(0,0,2,0), 
         col = high_contrast_pal)

dev.off()

cat(sprintf("[✓] Pipeline complete! Outputs saved to ./%s/\n", OUT_DIR))
