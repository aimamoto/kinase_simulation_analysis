# ==============================================================================
# Phase 8: Differential Allosteric Driver Analysis (Dual-Mode)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("[!] Required arguments: OUT_DIR TARGET_A TARGET_B GROUP_COLUMN")
}

OUT_DIR    <- args[1]
# Safely decode the newlines passed from Python
TARGET_A   <- gsub("___", "\n", args[2])
TARGET_B   <- gsub("___", "\n", args[3])
GROUP_COL  <- args[4]  # e.g., "Macro_State", "Condition", or "Group_Name"

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(ggrepel)
})

# ==============================================================================
# --- CONFIG: PLOT DIMENSIONS & AESTHETICS ---
# ==============================================================================
PDF_WIDTH             <- 10.0   # Width of the output PDF in inches
PDF_HEIGHT            <- 8.0    # Height of the output PDF in inches

POINT_SIZE            <- 3.0    # Size of the dots in the volcano plot
LABEL_SIZE            <- 4.5    # Size of the text labels (Note: in mm, not pts. 4.5 ≈ 12pt)

TITLE_FONT_SIZE       <- 14     # Font size for the main plot title
SUBTITLE_FONT_SIZE    <- 11     # Font size for the subtitle
AXIS_TITLE_FONT_SIZE  <- 12     # Font size for X and Y axis labels
AXIS_TEXT_FONT_SIZE   <- 10     # Font size for the tick numbers
LEGEND_FONT_SIZE      <- 11     # Font size for the legend text
# ==============================================================================

# Load the master metadata generated in Phase 7
master_csv <- file.path(OUT_DIR, "Phase7_Complete_Structural_Metadata.csv")
if (!file.exists(master_csv)) {
  stop(sprintf("[!] Cannot find %s", master_csv))
}

df <- read_csv(master_csv, show_col_types = FALSE)

# Check if the requested grouping column exists
if (!(GROUP_COL %in% colnames(df))) {
  stop(sprintf("[!] Column '%s' not found in metadata.", GROUP_COL))
}

# Filter exactly to the two targets requested
sub_df <- df %>% filter(!!sym(GROUP_COL) %in% c(TARGET_A, TARGET_B))

if (nrow(sub_df) < 10) {
  stop(sprintf("[!] Insufficient data to compare %s and %s.", TARGET_A, TARGET_B))
}

# Identify all valid distance metrics dynamically
dist_cols <- grep("_Dist$", colnames(sub_df), value = TRUE)
sub_df <- sub_df %>% mutate(across(all_of(dist_cols), as.numeric))

# Keep only columns that have enough data in BOTH groups to run stats
valid_cols <- c()
for(col in dist_cols) {
  counts <- sub_df %>% drop_na(all_of(col)) %>% count(!!sym(GROUP_COL))
  if (nrow(counts) == 2 && all(counts$n >= 5)) {
    valid_cols <- c(valid_cols, col)
  }
}

if (length(valid_cols) == 0) {
  stop("[!] No valid metrics found with sufficient N in both groups.")
}

# --- STATISTICAL TESTING ---
# Convert to long format for rstatix
df_long <- sub_df %>%
  select(Simulation_ID, !!sym(GROUP_COL), all_of(valid_cols)) %>%
  pivot_longer(cols = all_of(valid_cols), names_to = "Metric", values_to = "Distance") %>%
  drop_na()

# 1. Wilcoxon Rank Sum (FDR corrected)
stats_res <- df_long %>%
  group_by(Metric) %>%
  wilcox_test(as.formula(paste("Distance ~", GROUP_COL))) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  select(Metric, p, p.adj, p.adj.signif)

# 2. Cohen's d (Effect Size)
eff_res <- df_long %>%
  group_by(Metric) %>%
  cohens_d(as.formula(paste("Distance ~", GROUP_COL))) %>%
  select(Metric, effsize, magnitude)

# 3. Directionality (Log2 Fold Change equivalent for distances)
means <- df_long %>%
  group_by(Metric, !!sym(GROUP_COL)) %>%
  summarise(Mean_Dist = mean(Distance), .groups = "drop") %>%
  pivot_wider(names_from = !!sym(GROUP_COL), values_from = Mean_Dist)

# Merge everything into a Volcano DataFrame
volcano_df <- stats_res %>%
  left_join(eff_res, by = "Metric") %>%
  left_join(means, by = "Metric") %>%
  mutate(
    Log10_P = -log10(p.adj),
    # To determine direction, we subtract Target A from Target B
    # If Diff is positive, Target B is wider/expanded. If negative, Target B is compressed.
    Diff = .data[[TARGET_B]] - .data[[TARGET_A]],
    Direction = if_else(Diff > 0, "Expanded", "Compressed"),
    # Convert absolute effect size to a signed effect size for the X-axis
    Signed_EffSize = if_else(Direction == "Expanded", effsize, -effsize)
  ) %>%
  mutate(
    Significance = case_when(
      p.adj < 0.05 & Signed_EffSize > 0.5 ~ "Significantly Expanded",
      p.adj < 0.05 & Signed_EffSize < -0.5 ~ "Significantly Compressed",
      TRUE ~ "Not Significant"
    )
  )

# --- CREATE VOLCANO PLOT ---
# Strip newlines and @ symbols for a safe OS filename
clean_a <- gsub("[^A-Za-z0-9_\\-]", "_", TARGET_A)
clean_b <- gsub("[^A-Za-z0-9_\\-]", "_", TARGET_B)

# Handle string wrapping for beautiful plotting titles
title_a <- str_replace_all(TARGET_A, "\n|@", " & ")
title_b <- str_replace_all(TARGET_B, "\n|@", " & ")

p_volcano <- ggplot(volcano_df, aes(x = Signed_EffSize, y = Log10_P, color = Significance)) +
  geom_point(alpha = 0.8, size = POINT_SIZE) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Significantly Compressed" = "#2171b5", 
                                "Significantly Expanded" = "#cb181d", 
                                "Not Significant" = "grey70")) +
  geom_text_repel(data = filter(volcano_df, Significance != "Not Significant"),
                  aes(label = str_replace(Metric, "_Dist", "")),
                  size = LABEL_SIZE, max.overlaps = 15, box.padding = 0.5) +
  labs(title = sprintf("Differential Allosteric Drivers:\n%s  vs  %s", title_a, title_b),
       subtitle = sprintf("Positive Effect Size = Expanded in %s", title_b),
       x = "Cohen's d (Effect Size)",
       y = "-Log10(FDR Adjusted p-value)") +
  theme_bw() + 
  theme(
    legend.position = "bottom", 
    plot.title = element_text(face = "bold", hjust = 0.5, size = TITLE_FONT_SIZE),
    plot.subtitle = element_text(hjust = 0.5, size = SUBTITLE_FONT_SIZE),
    axis.title = element_text(size = AXIS_TITLE_FONT_SIZE),
    axis.text = element_text(size = AXIS_TEXT_FONT_SIZE),
    legend.title = element_text(size = LEGEND_FONT_SIZE),
    legend.text = element_text(size = LEGEND_FONT_SIZE)
  )

# Save outputs
dir.create(file.path(OUT_DIR, "Phase8_Volcanos"), showWarnings = FALSE)
ggsave(file.path(OUT_DIR, "Phase8_Volcanos", sprintf("Volcano_%s_vs_%s.pdf", clean_a, clean_b)), 
       p_volcano, width = PDF_WIDTH, height = PDF_HEIGHT)
write_csv(volcano_df, file.path(OUT_DIR, "Phase8_Volcanos", sprintf("Stats_%s_vs_%s.csv", clean_a, clean_b)))

cat("Done!\n")
