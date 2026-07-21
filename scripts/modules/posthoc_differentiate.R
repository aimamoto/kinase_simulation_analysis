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
GROUP_COL  <- args[4]  
EFF_METRIC <- if (length(args) > 4) args[5] else "wilcox" # "cohens_d" or "wilcox"

# --- Dependencies & Setup Safety Net ---
required_packages <- c("tidyverse", "rstatix", "ggrepel")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(ggrepel)
})

# ==============================================================================
# --- CONFIG: PLOT DIMENSIONS & AESTHETICS ---
# ==============================================================================
CFG <- list(
  pdf_width = 10.0,
  pdf_height = 8.0,
  point_size = 3.0,
  label_size = 4.5,
  title_font = 14,
  subtitle_font = 11,
  axis_title = 12,
  axis_text = 10,
  legend_font = 11
)
# ==============================================================================

master_csv <- file.path(OUT_DIR, "Phase7_Complete_Structural_Metadata.csv")
if (!file.exists(master_csv)) stop(sprintf("[!] Cannot find %s", master_csv))

df <- read_csv(master_csv, show_col_types = FALSE)

if (!(GROUP_COL %in% colnames(df))) stop(sprintf("[!] Column '%s' not found in metadata.", GROUP_COL))

sub_df <- df %>% filter(!!sym(GROUP_COL) %in% c(TARGET_A, TARGET_B))
if (nrow(sub_df) < 10) stop(sprintf("[!] Insufficient data to compare %s and %s.", TARGET_A, TARGET_B))

dist_cols <- grep("_Dist$", colnames(sub_df), value = TRUE)
sub_df <- sub_df %>% mutate(across(all_of(dist_cols), as.numeric))

valid_cols <- c()
for(col in dist_cols) {
  counts <- sub_df %>% drop_na(all_of(col)) %>% count(!!sym(GROUP_COL))
  # Require both groups (n>=5) AND non-zero variance across the subset. A column
  # that is constant everywhere carries no differential signal and makes
  # stats::wilcox.test abort ("missing value where TRUE/FALSE needed"); a column
  # constant *within* a group but differing *between* groups still has non-zero
  # overall variance, so it is kept (that is the perfect-discriminator case).
  col_var <- var(sub_df[[col]], na.rm = TRUE)
  if (nrow(counts) == 2 && all(counts$n >= 5) && is.finite(col_var) && col_var > 0)
    valid_cols <- c(valid_cols, col)
}
if (length(valid_cols) == 0) stop("[!] No valid metrics found with sufficient N in both groups.")

# --- STATISTICAL TESTING ---
# Force the grouping variable into a strict Factor to prevent ordering bugs
df_long <- sub_df %>%
  select(Simulation_ID, !!sym(GROUP_COL), all_of(valid_cols)) %>%
  pivot_longer(cols = all_of(valid_cols), names_to = "Metric", values_to = "Distance") %>%
  drop_na() %>%
  mutate(!!sym(GROUP_COL) := factor(!!sym(GROUP_COL), levels = c(TARGET_A, TARGET_B)))

# We keep statistic (W), n1, and n2 explicitly so we can do bulletproof math.
# stats::wilcox.test can abort ("missing value where TRUE/FALSE needed") on
# degenerate rank distributions (within-group-constant / perfect-separation with
# ties). Test each metric independently and, only if the default (exact) path
# errors, fall back to the normal approximation (exact = FALSE) so one
# pathological column cannot kill the whole comparison. Well-behaved columns keep
# their default p-values unchanged; a metric that fails both attempts is dropped
# (its p is uncomputable) rather than crashing the run.
wilcox_formula <- as.formula(paste("Distance ~", GROUP_COL))

# Last-resort manual Mann-Whitney U with a tie-corrected normal approximation.
# Handles degenerate rank distributions (e.g. a perfect discriminator: constant
# within each group but different between them) that make rstatix::wilcox_test
# abort even with exact = FALSE. Returns the same W (=U1)/n1/n2/p columns so the
# downstream rank-biserial math and BH adjustment are unaffected, and such
# columns are reported rather than silently dropped.
manual_wilcox <- function(dat) {
  a <- dat$Distance[dat[[GROUP_COL]] == TARGET_A]
  b <- dat$Distance[dat[[GROUP_COL]] == TARGET_B]
  n1 <- length(a); n2 <- length(b); N <- n1 + n2
  if (n1 < 1 || n2 < 1) return(NULL)
  rk  <- rank(c(a, b))
  U1  <- sum(rk[seq_len(n1)]) - n1 * (n1 + 1) / 2
  tab <- table(c(a, b)); tie <- sum(tab^3 - tab)
  sigma2 <- (n1 * n2 / 12) * ((N + 1) - tie / (N * (N - 1)))
  if (!is.finite(sigma2) || sigma2 <= 0) return(NULL)   # no rank variance at all
  z0 <- U1 - n1 * n2 / 2
  z  <- (z0 - sign(z0) * 0.5) / sqrt(sigma2)   # continuity correction, matching stats::wilcox.test
  tibble(.y. = "Distance", group1 = TARGET_A, group2 = TARGET_B,
         n1 = n1, n2 = n2, statistic = U1, p = 2 * pnorm(-abs(z)))
}

safe_wilcox <- function(dat) {
  tryCatch(
    wilcox_test(dat, wilcox_formula),
    error = function(e) tryCatch(
      wilcox_test(dat, wilcox_formula, exact = FALSE),
      error = function(e2) manual_wilcox(dat)
    )
  )
}

stats_res <- unique(df_long$Metric) %>%
  purrr::map(function(mn) {
    res <- safe_wilcox(dplyr::filter(df_long, Metric == mn))
    if (!is.null(res)) res$Metric <- mn
    res
  }) %>%
  purrr::compact() %>%
  dplyr::bind_rows()

if (nrow(stats_res) == 0) stop("[!] No metrics could be statistically tested for this comparison.")

stats_res <- stats_res %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

if (EFF_METRIC == "cohens_d") {
  # Manual, bulletproof Cohen's d computed directly from group moments, using the
  # SAME denominator as rstatix::cohens_d -- sqrt(mean(var1, var2)) -- so values
  # are numerically identical to the previous implementation on well-behaved
  # columns. rstatix::cohens_d errors internally (as_tidy_stat/tidy) on
  # degenerate zero-variance columns; computing it by hand keeps the pipeline
  # alive and yields clean Inf (perfect separation) / NaN (no difference) that
  # the finite-only axis bounds + clamped plot copy below handle gracefully.
  eff_res <- df_long %>%
    group_by(Metric, !!sym(GROUP_COL)) %>%
    summarise(.m = mean(Distance), .v = var(Distance), .groups = "drop_last") %>%
    summarise(
      m1 = dplyr::first(.m), m2 = dplyr::last(.m),
      v1 = dplyr::first(.v), v2 = dplyr::last(.v),
      .groups = "drop"
    ) %>%
    mutate(effsize = abs((m1 - m2) / sqrt((v1 + v2) / 2))) %>%
    select(Metric, effsize)
  x_axis_label <- "Cohen's d (Effect Size)"
  eff_threshold <- 0.5
} else {
  # INDESTRUCTIBLE NON-PARAMETRIC EFFECT SIZE (Rank-Biserial Correlation)
  # derived directly from the Mann-Whitney U (statistic) to bypass rstatix/coin tie-crashes
  eff_res <- stats_res %>%
    mutate(effsize = abs((2 * statistic) / (n1 * n2) - 1)) %>%
    select(Metric, effsize)
  x_axis_label <- "Rank-Biserial Correlation (r)"
  eff_threshold <- 0.3
}

# Simplify the stats df for joining
stats_clean <- stats_res %>% select(Metric, p, p.adj, p.adj.signif)

means <- df_long %>%
  group_by(Metric, !!sym(GROUP_COL)) %>%
  summarise(Mean_Dist = mean(Distance), .groups = "drop") %>%
  pivot_wider(names_from = !!sym(GROUP_COL), values_from = Mean_Dist)

volcano_df <- stats_clean %>%
  left_join(eff_res, by = "Metric") %>%
  left_join(means, by = "Metric") %>%
  mutate(
    Log10_P = -log10(p.adj),
    Diff = .data[[TARGET_B]] - .data[[TARGET_A]],
    Direction = if_else(Diff > 0, "Expanded", "Compressed"),
    Signed_EffSize = if_else(Direction == "Expanded", effsize, -effsize),
    Significance = case_when(
      p.adj < 0.05 & Signed_EffSize > eff_threshold ~ "Significantly Expanded",
      p.adj < 0.05 & Signed_EffSize < -eff_threshold ~ "Significantly Compressed",
      TRUE ~ "Not Significant"
    )
  )

# --- CREATE VOLCANO PLOT ---
clean_a <- gsub("[^A-Za-z0-9_\\-]", "_", TARGET_A)
clean_b <- gsub("[^A-Za-z0-9_\\-]", "_", TARGET_B)

title_a <- str_replace_all(TARGET_A, "\n|@", " & ")
title_b <- str_replace_all(TARGET_B, "\n|@", " & ")

# Force Symmetric X-Axis for Visual Balance.
# Cohen's d is Inf on zero-variance columns (pooled SD = 0), and p.adj can
# underflow to 0 (=> Log10_P = Inf). na.rm drops NA/NaN but NOT Inf, so derive
# axis bounds from finite values only (never feed infinite limits to ggsave),
# then plot a clamped copy so such points render at the margin instead of
# crashing or silently vanishing. volcano_df (the CSV) keeps the true values.
finite_eff <- volcano_df$Signed_EffSize[is.finite(volcano_df$Signed_EffSize)]
max_x <- if (length(finite_eff) > 0) max(abs(finite_eff)) * 1.1 else 1.0
if (!is.finite(max_x) || max_x == 0) max_x <- 1.0

finite_p <- volcano_df$Log10_P[is.finite(volcano_df$Log10_P)]
max_y <- if (length(finite_p) > 0) max(finite_p) * 1.1 else 1.0

plot_df <- volcano_df %>%
  mutate(
    Signed_EffSize = pmin(pmax(Signed_EffSize, -max_x), max_x),
    Log10_P        = pmin(Log10_P, max_y)
  )

p_volcano <- ggplot(plot_df, aes(x = Signed_EffSize, y = Log10_P, color = Significance)) +
  geom_point(alpha = 0.8, size = CFG$point_size) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-eff_threshold, eff_threshold), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Significantly Compressed" = "#2171b5", 
                                "Significantly Expanded" = "#cb181d", 
                                "Not Significant" = "grey70")) +
  geom_text_repel(data = filter(plot_df, Significance != "Not Significant"),
                  aes(label = str_replace(Metric, "_Dist", "")),
                  size = CFG$label_size, max.overlaps = 15, box.padding = 0.5) +
  scale_x_continuous(limits = c(-max_x, max_x)) +
  labs(title = sprintf("Differential Allosteric Drivers:\n%s  vs  %s", title_a, title_b),
       subtitle = sprintf("Positive Effect Size = Expanded in %s", title_b),
       x = x_axis_label,
       y = "-Log10(FDR Adjusted p-value)") +
  theme_bw() + 
  theme(
    legend.position = "bottom", 
    plot.title = element_text(face = "bold", hjust = 0.5, size = CFG$title_font),
    plot.subtitle = element_text(hjust = 0.5, size = CFG$subtitle_font),
    axis.title = element_text(size = CFG$axis_title),
    axis.text = element_text(size = CFG$axis_text),
    legend.title = element_text(size = CFG$legend_font),
    legend.text = element_text(size = CFG$legend_font)
  )

dir.create(file.path(OUT_DIR, "Phase8_Volcanos"), showWarnings = FALSE)
ggsave(file.path(OUT_DIR, "Phase8_Volcanos", sprintf("Volcano_%s_vs_%s.pdf", clean_a, clean_b)), 
       p_volcano, width = CFG$pdf_width, height = CFG$pdf_height)
write_csv(volcano_df, file.path(OUT_DIR, "Phase8_Volcanos", sprintf("Stats_%s_vs_%s.csv", clean_a, clean_b)))

cat("Done!\n")