# ==============================================================================
# Kinase Structural State Analysis - Unified Clustering Engine
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
TARGET_TYPE <- if (length(args) > 0) toupper(args[1]) else stop("Target Type required.")
CLUSTER_METHOD <- if (length(args) > 1) tolower(args[2]) else "gmm"
FILE_LIST_PATH <- if (length(args) > 2) args[3] else stop("CSV list file required.")

if (!(CLUSTER_METHOD %in% c("kmeans", "gmm"))) {
  stop("[!] Invalid CLUSTER_METHOD. Must be 'kmeans' or 'gmm'.")
}

METADATA_FILE <- "experiment.csv" 
HEATMAP_COLS <- 4

# Read dynamic output dir passed from Python wrapper
custom_out <- Sys.getenv("CUSTOM_OUT_DIR")
OUT_DIR <- if(nzchar(custom_out)) custom_out else sprintf("plots_and_stats_%s_%s", TARGET_TYPE, toupper(CLUSTER_METHOD))

TARGET_CATEGORIES <- list(
  State = "Active (BLAminus)", C_Helix = "In", R_Spine = "Intact", Spatial = "DFGin"
)

# --- Dependencies & Setup ---
required_packages <- c("tidyverse", "rstatix", "patchwork", "corrplot", "factoextra", 
                       "RColorBrewer", "cluster", "mclust", "plotly", "htmlwidgets", "pals")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(patchwork)
  library(corrplot)
  library(factoextra)
  library(RColorBrewer)
  library(pals)
  library(ggrepel)
})

# --- HELPER: TEXT WRAPPING FOR LONG PTM LABELS ---
wrap_text <- function(x, w = 42) {
  unname(sapply(x, function(s) {
    lines <- unlist(strsplit(as.character(s), "\n"))
    paste(sapply(lines, function(l) paste(strwrap(l, width = w), collapse = "\n")), collapse = "\n")
  }))
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[*] Pipeline Initialized.\n    Target: %s\n    Method: %s\n    Output: ./%s/\n", 
            TARGET_TYPE, toupper(CLUSTER_METHOD), OUT_DIR))

PLOT_THEME <- theme_classic() + 
  theme(text = element_text(size = 14, color = "black"), axis.text = element_text(color = "black"),
        legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))

# ==============================================================================
# DATA LOADING & METADATA MAPPING (SYSTEM-LEVEL)
# ==============================================================================
if(!file.exists(FILE_LIST_PATH)) stop(sprintf("[!] Cannot find file list at %s", FILE_LIST_PATH))
all_csvs <- readLines(FILE_LIST_PATH)
all_csvs <- all_csvs[all_csvs != ""]
if(length(all_csvs) == 0) stop("[!] No CSV files provided in the list.")

df <- map_dfr(all_csvs, ~read_csv(.x, col_types = cols(.default = "c"), 
                                  na = c("", "NA", "N/A", "None"), 
                                  show_col_types = FALSE, 
                                  name_repair = "unique_quiet")) %>%
  select(-matches("\\.\\.\\.\\d+$")) %>%
  distinct(Simulation_ID, Chain, .keep_all = TRUE)

# ------------------------------------------------------------------------------
# 1. ONTOLOGICAL METADATA MAPPING OR AUTO-DISCOVERY
# ------------------------------------------------------------------------------
# Helper: Strips AF3 seed/sample tags so 100 directories cleanly collapse to 1 group label
clean_sim_name <- function(name) {
  if(is.na(name)) return(NA_character_)
  s <- gsub("_?(seed|sample|model|rank|pred|unrelaxed|relaxed)[-_]?[0-9]*.*$", "", name, ignore.case = TRUE)
  s <- gsub("_?(unrelaxed|relaxed)$", "", s, ignore.case = TRUE)
  return(s)
}

# Helper: Auto-formats a raw directory string into a clean plot label WITHOUT a-/b- prefixes
format_auto_label <- function(sim_id) {
  parts <- unlist(strsplit(sim_id, "_"))
  if(length(parts) == 1) return(sub("^[abAB]-", "", sim_id))
  
  cond <- parts[length(parts)]
  prots <- parts[1:(length(parts)-1)]
  prots <- sub("^[abAB]-", "", prots) # Strip prefixes cosmetically
  paste0(paste(prots, collapse="\n"), "\n[", cond, "]")
}

df <- df %>% mutate(Clean_Sim_ID = sapply(Simulation_ID, clean_sim_name))

if (file.exists(METADATA_FILE)) {
  cat(sprintf("\n[*] Found '%s'. Applying Strict Ontological Mapping (Hyphen/Underscore logic)...\n", METADATA_FILE))
  meta <- read_csv(METADATA_FILE, show_col_types = FALSE)
  colnames(meta) <- tolower(trimws(colnames(meta)))
  
  expected_cols <- c("chain_a", "mutation_a", "ptm_a", "chain_b", "mutation_b", "ptm_b")
  for (col in expected_cols) if (!(col %in% colnames(meta))) meta[[col]] <- ""
  
  build_chain <- function(chain, mut, ptm) {
    base <- c(chain, mut, ptm)
    base <- base[base != "" & !is.na(base)]
    paste(base, collapse="-")
  }
  
  cond_col <- if("condition" %in% colnames(meta)) "condition" else "condition_a"
  
  meta <- meta %>%
    rowwise() %>%
    mutate(
      clean_a = tolower(build_chain(chain_a, mutation_a, ptm_a)),
      clean_b = tolower(build_chain(chain_b, mutation_b, ptm_b)),
      cond_val = tolower(get(cond_col)),
      
      # Mathematically generate the expected Simulation ID based on user convention
      expected_id_1 = if_else(clean_b == "", paste0(clean_a, "_", cond_val), paste0("a-", clean_a, "_b-", clean_b, "_", cond_val)),
      expected_id_2 = if_else(clean_b == "", paste0("a-", clean_a, "_", cond_val), paste0("b-", clean_b, "_a-", clean_a, "_", cond_val)),
      
      # Generate a perfectly structured multi-line plot label WITHOUT a-/b- prefixes
      Condition_Label = if_else(clean_b == "", 
                                paste0(clean_a, "\n[", cond_val, "]"), 
                                paste0(clean_a, "\n", clean_b, "\n[", cond_val, "]"))
    ) %>% ungroup()
    
  # Logic Engine: Maps AF3 directory -> experiment.csv group
  match_metadata <- function(sim_id, meta_df) {
    sim_lower <- tolower(sim_id)
    
    # 1. Direct Convention Check (Perfect match)
    for(i in 1:nrow(meta_df)) {
      if(sim_lower == meta_df$expected_id_1[i] || sim_lower == meta_df$expected_id_2[i]) {
        return(meta_df$Condition_Label[i])
      }
    }
    
    # 2. Substring Convention Check (e.g., if trailing characters were missed)
    for(i in 1:nrow(meta_df)) {
      if(startsWith(sim_lower, meta_df$expected_id_1[i]) || startsWith(sim_lower, meta_df$expected_id_2[i])) {
        return(meta_df$Condition_Label[i])
      }
    }
    
    # 3. Robust Logical Token Match (Fallback)
    best_match <- format_auto_label(sim_id) # Default to clean auto-label if no match
    max_score <- -1
    for(i in 1:nrow(meta_df)) {
      row <- meta_df[i, ]
      req_tokens <- c(unlist(strsplit(row$clean_a, "-")), unlist(strsplit(row$clean_b, "-")), unlist(strsplit(row$cond_val, "-")))
      req_tokens <- req_tokens[req_tokens != ""]
      
      if(length(req_tokens) > 0 && all(sapply(req_tokens, function(x) grepl(x, sim_lower, fixed=TRUE)))) {
         if(length(req_tokens) > max_score) {
           max_score <- length(req_tokens)
           best_match <- row$Condition_Label
         }
      }
    }
    return(best_match)
  }
  
  df <- df %>% mutate(Condition_designed = sapply(Clean_Sim_ID, function(x) match_metadata(x, meta)))
  df <- df %>% filter(Condition_designed %in% meta$Condition_Label)
  
} else {
  cat(sprintf("\n[*] No '%s' found. Activating Auto-Discovery Mode...\n", METADATA_FILE))
  df <- df %>% mutate(Condition_designed = sapply(Clean_Sim_ID, format_auto_label))
}

if(nrow(df) == 0) stop("[!] Metadata filtering removed all data. Ensure your directory strings match experiment.csv components.")

# ==============================================================================
# Phase 1: SYSTEM-LEVEL LIGAND OCCUPANCY (TUG-OF-WAR)
# ==============================================================================
if (any(grepl("1atp", df$Condition_designed, ignore.case = TRUE))) {
  cat("\n[*] Executing Phase 1: System-Level Ligand Occupancy (Tug-of-War)...\n")
  
  occ_df <- df %>%
    group_by(Simulation_ID, Condition_designed) %>%
    summarise(
      Bound_Types = paste(sort(Type[!grepl("no ligand", C_Spine, ignore.case=TRUE)]), collapse = " + "),
      .groups = "drop"
    ) %>%
    mutate(Bound_State = if_else(Bound_Types == "", "Apo (Neither Bound)", paste(Bound_Types, "Bound")))
  
  occ_summary <- occ_df %>%
    mutate(Condition_clean = gsub("\\s*\\[.*?\\]", "", Condition_designed)) %>%
    group_by(Condition_clean, Bound_State) %>%
    summarise(Count = n(), .groups="drop") %>%
    group_by(Condition_clean) %>%
    mutate(Percent = Count / sum(Count) * 100)
  
  p_occ <- ggplot(occ_summary, aes(x = Condition_clean, y = Percent, fill = Bound_State)) +
    geom_bar(stat="identity", color="black", width=0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Dimer-Level Ligand Occupancy (AF3 Physical Coordinates)",
         subtitle = "Illustrates 1ATP Tug-of-war preference and validates Apo/Holo controls",
         x = "Condition", y = "Proportion of Models (%)") +
    PLOT_THEME + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10),
                       legend.title = element_blank())
  
  ggsave(file.path(OUT_DIR, "Phase1_Ligand_Occupancy.pdf"), p_occ, width=10, height=7)
  write_csv(occ_summary, file.path(OUT_DIR, "Phase1_Ligand_Occupancy_Stats.csv"))
  cat("    -> [✓] Phase 1 Complete. Check 'Phase1_Ligand_Occupancy.pdf'.\n")
} else {
  cat("\n[*] Phase 1 Skipped: No '1atp' (competitive) conditions detected in the filtered cohort.\n")
}

# ==============================================================================
# TARGET ISOLATION
# ==============================================================================
cat(sprintf("\n[*] Isolating %s Target Data for Structural Analysis...\n", TARGET_TYPE))

# Safely tolerate A- or B- prefixes natively assigned by extract_fasta.py
df <- df %>% filter(grepl(paste0("^(?:[A-Z]-)?", TARGET_TYPE), toupper(Type)))

df <- df %>% mutate(
  Actual_Target_Apo = grepl("no ligand", C_Spine, ignore.case = TRUE),
  Physical_State = if_else(Actual_Target_Apo, "Apo", "Holo")
)

# Strips the dynamically generated [condition] block to append the physical state
df$Condition_clean <- gsub("\\s*\\[.*?\\]", "", df$Condition_designed)
df <- df %>% mutate(
  Condition_physical = paste0(Condition_clean, "\n[", TARGET_TYPE, " = ", Physical_State, "]")
)

df$Condition <- factor(df$Condition_physical)
group_counts <- df %>% group_by(Condition) %>% summarise(n = n(), .groups="drop")
valid_groups <- group_counts %>% filter(n >= 5) %>% pull(Condition)

if(length(valid_groups) == 0) stop("[!] FATAL: No condition groups have enough valid simulations (n >= 5) to run statistics.")

df <- df %>% filter(Condition %in% valid_groups) %>% droplevels()

unique_conditions <- levels(df$Condition)
num_conditions <- length(unique_conditions)
safe_colors <- brewer.pal(max(3, 9), "Set1")[-6] 
if (num_conditions <= 8) {
  CUSTOM_COLORS <- setNames(safe_colors[1:num_conditions], unique_conditions)
} else {
  CUSTOM_COLORS <- setNames(colorRampPalette(safe_colors)(num_conditions), unique_conditions)
}

all_dist_cols <- grep("_Dist$|^RAF_", colnames(df), value = TRUE)
target_cols <- c("Phi_D", "Psi_D", all_dist_cols)
df <- df %>% mutate(across(any_of(target_cols), ~ suppressWarnings(as.numeric(as.character(.)))))

ligand_cols <- grep("ATP|Mg", all_dist_cols, ignore.case = TRUE, value = TRUE)
if(length(ligand_cols) > 0) {
  df <- df %>% mutate(across(any_of(ligand_cols), ~ if_else(Actual_Target_Apo, NA_real_, .)))
}

viable_dist_summary <- df %>% select(Condition, all_of(all_dist_cols)) %>%
  group_by(Condition) %>% summarise(across(everything(), ~sum(!is.na(.)) >= 5)) %>% select(-Condition)

universal_dist_cols <- names(which(sapply(viable_dist_summary, all)))
any_dist_cols <- names(which(sapply(viable_dist_summary, any)))
multi_dist_cols <- names(which(sapply(viable_dist_summary, sum) >= 2))

if(length(universal_dist_cols) < 2) {
  stop("[!] FATAL: Less than 2 universal distance metrics found. PCA cannot be computed.")
}

# ==============================================================================
# Phase 2: Macro-State Conformational Shifts
# ==============================================================================
cat("[*] Executing Phase 2: Categorical Macro-States...\n")

plot_categorical <- function(v_name, title) {
  df %>% group_by(Condition, !!sym(v_name)) %>% summarise(n = n(), .groups="drop") %>% 
    group_by(Condition) %>% mutate(Percent = n / sum(n) * 100) %>%
    ggplot(aes(x = Condition, y = Percent, fill = !!sym(v_name))) +
    geom_bar(stat = "identity", color = "black", width = 0.6) + scale_fill_brewer(palette = "Paired") +
    labs(title = sprintf("%s: %s", TARGET_TYPE, title), x = NULL, y = "Proportion (%)", fill = "State") +
    PLOT_THEME + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
}
p_phase2 <- (plot_categorical("State", "Global Conformation") | plot_categorical("C_Helix", "aC-Helix State")) / 
  (plot_categorical("R_Spine", "R-Spine Integrity") | plot_categorical("Spatial", "DFG Spatial State"))
ggsave(file.path(OUT_DIR, "Phase2_Macro_States.pdf"), p_phase2, width = 16, height = 12)

pairwise_fisher_results <- tibble(Structural_Feature=character(), Group_1=character(), Group_2=character(), Target_State=character(), p_value=numeric())
condition_pairs <- combn(unique_conditions, 2, simplify = FALSE)
for(cat_var in names(TARGET_CATEGORIES)) {
  for(pair in condition_pairs) {
    sub_df <- df %>% filter(Condition %in% c(pair[1], pair[2])) %>%
      mutate(Binary_State = factor(if_else(!!sym(cat_var) == TARGET_CATEGORIES[[cat_var]], TARGET_CATEGORIES[[cat_var]], "Other")))
    tbl <- table(droplevels(sub_df$Condition), sub_df$Binary_State)
    if(nrow(tbl) == 2 && ncol(tbl) == 2 && sum(colSums(tbl) > 0) >= 2) {
      res <- tryCatch({ fisher.test(tbl) }, error = function(e) { fisher.test(tbl, simulate.p.value = TRUE, B = 2000) })
      pairwise_fisher_results <- pairwise_fisher_results %>% add_row(Structural_Feature=cat_var, Group_1=pair[1], Group_2=pair[2], Target_State=TARGET_CATEGORIES[[cat_var]], p_value=res$p.value)
    }
  }
}
write_csv(pairwise_fisher_results %>% group_by(Structural_Feature) %>% mutate(p.adj=p.adjust(p_value, "BH")) %>% add_significance("p.adj"), file.path(OUT_DIR, "Phase2_Pairwise_Categorical_Stats.csv"))

# ==============================================================================
# Phase 3: 2D Phase Space
# ==============================================================================
cat("[*] Executing Phase 3: 2D Phase Space...\n")
p_dihedral <- df %>% drop_na(Phi_D, Psi_D) %>% ggplot(aes(Phi_D, Psi_D, color=Condition, fill=Condition)) + geom_point(alpha=0.6) + geom_density_2d(alpha=0.8) + scale_color_manual(values=CUSTOM_COLORS, labels=function(x) wrap_text(x, 40)) + scale_fill_manual(values=CUSTOM_COLORS, labels=function(x) wrap_text(x, 40)) + labs(title="DFG Dihedral Phase Space", x=expression(Phi~"(°)"), y=expression(Psi~"(°)")) + PLOT_THEME
p_dunbrack <- df %>% drop_na(D1_Dist, D2_Dist) %>% ggplot(aes(D1_Dist, D2_Dist, color=Condition, fill=Condition)) + geom_point(alpha=0.6) + geom_density_2d(alpha=0.8) + scale_color_manual(values=CUSTOM_COLORS, labels=function(x) wrap_text(x, 40)) + scale_fill_manual(values=CUSTOM_COLORS, labels=function(x) wrap_text(x, 40)) + labs(title="Dunbrack Coordinates", x="D1 (Å)", y="D2 (Å)") + PLOT_THEME
suppressWarnings(ggsave(file.path(OUT_DIR, "Phase3_2D_PhaseSpace.pdf"), p_dihedral | p_dunbrack, width = 14, height = 6))

# ==============================================================================
# Phase 4: 1D Allosteric Stats
# ==============================================================================
cat("[*] Executing Phase 4: 1D Micro-Metrics & Stats...\n")
df_long <- df %>% select(Simulation_ID, Condition, all_of(any_dist_cols)) %>% pivot_longer(-c(Simulation_ID, Condition), names_to="Metric", values_to="Distance") %>% drop_na(Distance)

p_violins <- df_long %>% ggplot(aes(x=Condition, y=Distance, fill=Condition)) + 
  geom_violin(trim=FALSE, alpha=0.7, color=NA) + geom_boxplot(width=0.2, outlier.shape=NA, alpha=0.8) + 
  facet_wrap(~Metric, scales="free_y", ncol=4) + scale_fill_manual(values=CUSTOM_COLORS, labels=function(x) wrap_text(x, 40)) + 
  PLOT_THEME + theme(axis.text.x=element_blank())
ggsave(file.path(OUT_DIR, "Phase4_Allosteric_Distances.pdf"), p_violins, width = 16, height = 12)

df_stats <- df_long %>% filter(Metric %in% multi_dist_cols)
write_csv(df_stats %>% group_by(Metric) %>% kruskal_test(Distance ~ Condition) %>% adjust_pvalue(method = "BH") %>% add_significance("p.adj"), file.path(OUT_DIR, "Phase4_Kruskal_Summary.csv"))
write_csv(df_stats %>% group_by(Metric) %>% pairwise_wilcox_test(Distance ~ Condition, p.adjust.method="BH") %>% add_significance("p.adj"), file.path(OUT_DIR, "Phase4_Pairwise_Summary.csv"))

# ==============================================================================
# Phase 5: Network Coupling (PCA & Correlation)
# ==============================================================================
cat("[*] Executing Phase 5: Network Coupling (PCA & Correlation)...\n")
df_pca_complete <- df %>% drop_na(all_of(universal_dist_cols))

if(nrow(df_pca_complete) < 3) stop("[!] FATAL: Not enough complete data rows to perform PCA.")

pca_data <- df_pca_complete %>% select(all_of(universal_dist_cols))
pca_metadata <- df_pca_complete 

pca_res <- prcomp(pca_data, scale. = TRUE) 
suppressWarnings({
  p_pca <- fviz_pca_ind(pca_res, geom="point", col.ind=pca_metadata$Condition, palette=as.character(CUSTOM_COLORS), alpha=0.7, title=sprintf("%s: Shared Allosteric Variance", TARGET_TYPE)) + PLOT_THEME 
})
ggsave(file.path(OUT_DIR, "Phase5_PCA_Biplot.pdf"), p_pca, width = 10, height = 7)

high_contrast_pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)
n_cols <- min(HEATMAP_COLS, num_conditions)
n_rows <- ceiling(num_conditions / n_cols)

pdf(file.path(OUT_DIR, "Phase5_Correlation_Heatmaps.pdf"), width = 7 * n_cols, height = 7 * n_rows)
par(mfrow = c(n_rows, n_cols))

for (grp in unique_conditions) {
  df_grp <- df %>% filter(Condition == grp) %>% select(all_of(any_dist_cols)) %>% select(where(~sum(!is.na(.)) >= 3))
  if (nrow(df_grp) >= 3 && ncol(df_grp) >= 2) {
    cor_matrix <- suppressWarnings(cor(df_grp, use = "pairwise.complete.obs", method = "spearman"))
    cor_matrix[is.na(cor_matrix)] <- 0 
    
    title_str <- sprintf("%s Coupling\n(n=%d)", wrap_text(grp, 35), nrow(df_grp))
    
    corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black", tl.cex = 1.1, cl.cex = 1.1, addgrid.col = "white",
             title = title_str, mar = c(0,0,7,0), cex.main = 1.3, col = high_contrast_pal)
  } else {
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(0, 2), ylim = c(0, 2))
    text(1, 1.2, sprintf("%s Coupling", wrap_text(grp, 35)), font = 2, cex = 1.2)
    text(1, 0.8, sprintf("Insufficient Data\n(n=%d)", nrow(df_grp)), cex = 1.3, col = "#e74c3c")
  }
}
dev.off()

cat("    -> Quantifying Global Network Density & Differential Correlations...\n")

cor_matrices <- list()
n_counts <- list()
edge_weight_dists <- list() 
density_stats <- tibble(Condition = character(), N = numeric(), Global_Coupling_Score = numeric())

for (grp in unique_conditions) {
  df_grp <- df %>% filter(Condition == grp) %>% select(all_of(any_dist_cols)) %>% select(where(~sum(!is.na(.)) >= 3))
  if (nrow(df_grp) >= 3 && ncol(df_grp) >= 2) {
    c_mat <- suppressWarnings(cor(df_grp, use = "pairwise.complete.obs", method = "spearman"))
    c_mat[is.na(c_mat)] <- 0 
    cor_matrices[[grp]] <- c_mat
    n_counts[[grp]] <- nrow(df_grp)
    
    abs_edges <- abs(c_mat[upper.tri(c_mat)])
    edge_weight_dists[[grp]] <- abs_edges 
    
    mac <- mean(abs_edges)
    density_stats <- density_stats %>% add_row(Condition = grp, N = nrow(df_grp), Global_Coupling_Score = mac)
  }
}
write_csv(density_stats %>% arrange(desc(Global_Coupling_Score)), file.path(OUT_DIR, "Phase5_Global_Network_Density.csv"))

mac_stats <- tibble(Condition_1 = character(), Condition_2 = character(), 
                    MAC_1 = numeric(), MAC_2 = numeric(), p_value = numeric())

if(length(edge_weight_dists) >= 2) {
  mac_pairs <- combn(names(edge_weight_dists), 2, simplify = FALSE)
  for(mp in mac_pairs) {
    g1 <- mp[1]; g2 <- mp[2]
    edges1 <- edge_weight_dists[[g1]]; edges2 <- edge_weight_dists[[g2]]
    w_test <- tryCatch({ suppressWarnings(wilcox.test(edges1, edges2)) }, error = function(e) list(p.value=NA))
    mac_stats <- mac_stats %>% add_row(Condition_1 = g1, Condition_2 = g2, MAC_1 = mean(edges1), MAC_2 = mean(edges2), p_value = w_test$p.value)
  }
  mac_stats <- mac_stats %>% filter(!is.na(p_value))
  if(nrow(mac_stats) > 0) {
    mac_stats <- mac_stats %>% mutate(p.adj = p.adjust(p_value, "BH")) %>% add_significance("p.adj") %>% arrange(p.adj)
    write_csv(mac_stats, file.path(OUT_DIR, "Phase5_Global_Network_Density_Stats.csv"))
  }
}

diff_cor_stats <- tibble(Condition_1=character(), Condition_2=character(), Metric_1=character(), Metric_2=character(), 
                         R1=numeric(), R2=numeric(), Z_Score=numeric(), p_value=numeric())

if(length(cor_matrices) >= 2) {
  group_pairs <- combn(names(cor_matrices), 2, simplify = FALSE)
  for(pair in group_pairs) {
    g1 <- pair[1]; g2 <- pair[2]
    mat1 <- cor_matrices[[g1]]; mat2 <- cor_matrices[[g2]]
    n1 <- n_counts[[g1]]; n2 <- n_counts[[g2]]
    
    shared_metrics <- intersect(rownames(mat1), rownames(mat2))
    if(length(shared_metrics) >= 2 && n1 > 3 && n2 > 3) {
      metric_pairs <- combn(shared_metrics, 2, simplify = FALSE)
      for(mp in metric_pairs) {
        m_a <- mp[1]; m_b <- mp[2]
        r1 <- mat1[m_a, m_b]; r2 <- mat2[m_a, m_b]
        
        r1_safe <- min(max(r1, -0.999), 0.999); r2_safe <- min(max(r2, -0.999), 0.999)
        z1 <- 0.5 * log((1 + r1_safe) / (1 - r1_safe))
        z2 <- 0.5 * log((1 + r2_safe) / (1 - r2_safe))
        se_diff <- sqrt((1 / (n1 - 3)) + (1 / (n2 - 3)))
        z_score <- (z1 - z2) / se_diff
        p_val <- 2 * (1 - pnorm(abs(z_score)))
        
        diff_cor_stats <- diff_cor_stats %>% add_row(Condition_1=g1, Condition_2=g2, Metric_1=m_a, Metric_2=m_b, R1=r1, R2=r2, Z_Score=z_score, p_value=p_val)
      }
    }
  }
  
  if(nrow(diff_cor_stats) > 0) {
    diff_cor_stats <- diff_cor_stats %>% mutate(p.adj = p.adjust(p_value, "BH")) %>% add_significance("p.adj") %>% arrange(p.adj)
    write_csv(diff_cor_stats, file.path(OUT_DIR, "Phase5_Differential_Correlations.csv"))
    
    p_mac <- density_stats %>%
      ggplot(aes(x = reorder(Condition, Global_Coupling_Score), y = Global_Coupling_Score, fill = Condition)) +
      geom_col(color = "black", width = 0.7, alpha = 0.85) +
      coord_flip() + scale_fill_manual(values = CUSTOM_COLORS, labels = function(x) wrap_text(x, 40)) +
      scale_x_discrete(labels = function(x) wrap_text(x, 30)) +
      labs(title = sprintf("%s: Global Network Rigidity", TARGET_TYPE), subtitle = "Mean Absolute Correlation (MAC) of internal dynamics", x = NULL, y = "Global Coupling Score (MAC)") +
      PLOT_THEME + theme(legend.position = "none", axis.text.y = element_text(face="bold"))
    ggsave(file.path(OUT_DIR, "Phase5_Global_Network_Density.pdf"), p_mac, width = 10, height = 6)
    
    sig_pairs <- diff_cor_stats %>% filter(p.adj < 0.05) %>% distinct(Condition_1, Condition_2)
    if(nrow(sig_pairs) > 0) {
      pdf(file.path(OUT_DIR, "Phase5_Correlation_Shift_Volcanos.pdf"), width = 10, height = 8)
      for(i in 1:nrow(sig_pairs)) {
        c1 <- sig_pairs$Condition_1[i]; c2 <- sig_pairs$Condition_2[i]
        plot_data <- diff_cor_stats %>% filter(Condition_1 == c1 & Condition_2 == c2) %>%
          mutate(Delta_R = R2 - R1, Log10_P = -log10(p.adj),
                 Edge_Label = paste0(str_replace(Metric_1, "_Dist", ""), " : ", str_replace(Metric_2, "_Dist", "")),
                 Significance = case_when(p.adj < 0.05 & Delta_R > 0.4 ~ "Significantly Coupled", p.adj < 0.05 & Delta_R < -0.4 ~ "Significantly Decoupled", TRUE ~ "Not Significant"))
        
        title_1 <- wrap_text(str_replace_all(c1, "\n", " & "), 40)
        title_2 <- wrap_text(str_replace_all(c2, "\n", " & "), 40)
        
        p_volc <- ggplot(plot_data, aes(x = Delta_R, y = Log10_P, color = Significance)) +
          geom_point(alpha = 0.8, size = 3) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + geom_vline(xintercept = c(-0.4, 0.4), linetype = "dashed", color = "black") +
          scale_color_manual(values = c("Significantly Decoupled" = "#2171b5", "Significantly Coupled" = "#cb181d", "Not Significant" = "grey70")) +
          geom_text_repel(data = filter(plot_data, Significance != "Not Significant"), aes(label = Edge_Label), size = 3, max.overlaps = 15, box.padding = 0.5) +
          labs(title = "Network Rewiring", subtitle = sprintf("%s\n  vs  \n%s", title_1, title_2), x = expression(Delta*R~~(R[2] - R[1])), y = "-Log10(FDR Adjusted p-value)") +
          theme_bw() + theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust=0.5, size=14), plot.subtitle = element_text(hjust=0.5, size=11))
        print(p_volc)
      }
      dev.off()
    }
  }
}

# ==============================================================================
# Phase 6: DYNAMIC CLUSTERING ENGINE
# ==============================================================================
cat(sprintf("\n[*] Executing Phase 6: Unsupervised Discovery via %s...\n", toupper(CLUSTER_METHOD)))

get_distinct_colors <- function(k) {
  tab10_bold <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  tab10_light <- c("#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", 
                   "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
  tab20_smart <- c(tab10_bold, tab10_light)
  if (k <= 20) { return(tab20_smart[1:k]) } else { return(colorRampPalette(tab20_smart)(k)) }
}

if (CLUSTER_METHOD == "kmeans") {
  cluster_data <- pca_res$x[, 1:min(3, ncol(pca_res$x))]
  set.seed(42)
  suppressMessages(library(cluster))
  gap_stat <- clusGap(cluster_data, FUN = kmeans, nstart = 25, K.max = 8, B = 100)
  optimal_k <- max(2, maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method="globalSEmax"))
  km_final <- kmeans(cluster_data, centers = optimal_k, nstart = 50)
  pca_metadata$Macro_State <- factor(paste("State", km_final$cluster))
  
  suppressWarnings({
    p_gap <- fviz_gap_stat(gap_stat, maxSE=list(method="globalSEmax")) + PLOT_THEME
    ggsave(file.path(OUT_DIR, "Phase6_Statistical_Proof.pdf"), p_gap, width = 8, height = 6)
  })
} else if (CLUSTER_METHOD == "gmm") {
  cluster_data <- pca_res$x[, 1:min(4, ncol(pca_res$x))]
  set.seed(42)
  suppressPackageStartupMessages(library(mclust))
  gmm_res <- Mclust(cluster_data)
  optimal_k <- gmm_res$G
  pca_metadata$Macro_State <- factor(paste("State", gmm_res$classification))
  detach("package:mclust", unload=TRUE)
}

state_counts <- pca_metadata %>% group_by(Macro_State) %>% summarise(n=n(), .groups="drop")
threshold <- max(10, nrow(pca_metadata) * 0.05)
small_states <- state_counts %>% filter(n < threshold) %>% pull(Macro_State)

if(length(small_states) > 0 && length(unique(pca_metadata$Macro_State)) > 1) {
  centroids <- pca_metadata %>% group_by(Macro_State) %>% summarise(across(all_of(universal_dist_cols), mean)) %>% column_to_rownames("Macro_State")
  for (ss in small_states) {
    valid_targets <- setdiff(rownames(centroids), small_states)
    if(length(valid_targets) == 0) break
    dists <- as.matrix(dist(centroids))
    nearest <- names(which.min(dists[ss, valid_targets]))
    pca_metadata <- pca_metadata %>% mutate(Macro_State = as.character(Macro_State), Macro_State = ifelse(Macro_State == ss, nearest, Macro_State))
  }
  pca_metadata$Macro_State <- factor(pca_metadata$Macro_State)
  optimal_k <- length(unique(pca_metadata$Macro_State))
  cat(sprintf("    -> Consolidated micro-clusters. Final Meta-States: %d\n", optimal_k))
}

state_colors <- get_distinct_colors(optimal_k)
names(state_colors) <- levels(pca_metadata$Macro_State)

plot_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], State = pca_metadata$Macro_State)
hulls <- plot_df %>% group_by(State) %>% slice(chull(PC1, PC2)) %>% ungroup()

p_state_pca <- ggplot(plot_df, aes(x = PC1, y = PC2, color = State, fill = State)) +
  geom_point(alpha = 0.7, size = 2) + geom_polygon(data = hulls, alpha = 0.2, linewidth = 0) +
  scale_color_manual(values = state_colors) + scale_fill_manual(values = state_colors) +
  labs(title = sprintf("%s Meta-States", TARGET_TYPE)) + PLOT_THEME
ggsave(file.path(OUT_DIR, "Phase6_State_Clusters_PCA.pdf"), p_state_pca, width = 10, height = 7)

suppressPackageStartupMessages(library(plotly))
ax_style <- list(titlefont = list(size = 18, color = "black", family = "Arial"), tickfont = list(size = 14, color = "black", family = "Arial"))
p_3d <- plot_ly(data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2], PC3=pca_res$x[,3], State=pca_metadata$Macro_State), 
                x=~PC1, y=~PC2, z=~PC3, color=~State, colors=state_colors, type='scatter3d', mode='markers', marker=list(size=4, opacity=0.8, line=list(width=0))) %>%
  plotly::layout(title = list(text = sprintf("<b>%s Interactive 3D Phase Space</b>", toupper(CLUSTER_METHOD)), font = list(size = 24, color = "black")),
                 scene = list(xaxis = c(list(title = 'PC1'), ax_style), yaxis = c(list(title = 'PC2'), ax_style), zaxis = c(list(title = 'PC3'), ax_style)),
                 legend = list(title = list(text = "<b>Discovered States</b>", font = list(size = 22, color = "black")), font = list(size = 20, color = "black"), itemsizing = 'constant', itemwidth = 40), margin = list(t = 60))
suppressWarnings(htmlwidgets::saveWidget(p_3d, file.path(OUT_DIR, "Phase6_Interactive_3D_Space.html")))
detach("package:plotly", unload=TRUE)

p_state_comp <- pca_metadata %>% group_by(Macro_State, Condition) %>% summarise(n=n(), .groups="drop") %>% 
  ggplot(aes(x = Macro_State, y = n, fill = Condition)) + geom_bar(stat = "identity", color = "black", width = 0.6) + 
  scale_fill_manual(values = CUSTOM_COLORS, labels = function(x) wrap_text(x, 40)) + 
  labs(title = "Composition of Meta-Stable States", x = "Discovered State", y = "Number of Models") + PLOT_THEME
ggsave(file.path(OUT_DIR, "Phase6_State_Composition.pdf"), p_state_comp, width=10, height=7)

safe_cols <- c("Simulation_ID", "Condition_designed", "Condition_reviewed", "Macro_State")
pca_metadata_out <- pca_metadata %>% select(any_of(safe_cols))
write_csv(pca_metadata_out, file.path(OUT_DIR, "Phase6_State_Assignments.csv"))

# ==============================================================================
# Phase 7: Biological Signatures
# ==============================================================================
cat("\n[*] Executing Phase 7: Mapping States to Biological Features...\n")
for(cv in names(TARGET_CATEGORIES)) {
  pca_metadata[[paste0("Binary_", cv)]] <- factor(if_else(pca_metadata[[cv]]==TARGET_CATEGORIES[[cv]], TARGET_CATEGORIES[[cv]], "Other"), levels=c(TARGET_CATEGORIES[[cv]], "Other"))
}

plot_feat <- function(v_name, o_var) {
  pca_metadata %>% group_by(Macro_State, !!sym(v_name)) %>% summarise(n = n(), .groups="drop") %>% 
    group_by(Macro_State) %>% mutate(Percent=n/sum(n)*100) %>%
    ggplot(aes(x=Macro_State, y=Percent, fill=!!sym(v_name))) + geom_bar(stat="identity", color="black", width=0.6) + 
    scale_fill_manual(values=setNames(c("#2171b5", "#cccccc"), c(TARGET_CATEGORIES[[o_var]], "Other"))) +
    labs(title=sprintf("%s Identity", o_var), x=NULL, y="Proportion (%)", fill="Feature") + PLOT_THEME
}
ggsave(file.path(OUT_DIR, "Phase7_MacroState_Signatures.pdf"), (plot_feat("Binary_State","State")|plot_feat("Binary_C_Helix","C_Helix"))/(plot_feat("Binary_R_Spine","R_Spine")|plot_feat("Binary_Spatial","Spatial")), width=16, height=12)

macro_stats <- tibble(Structural_Feature=character(), Target_State=character(), p_value=numeric())
for(cv in names(TARGET_CATEGORIES)) {
  tbl_data <- pca_metadata %>% group_by(Macro_State, !!sym(paste0("Binary_", cv))) %>% summarise(Count = n(), .groups="drop") %>%
    pivot_wider(names_from = !!sym(paste0("Binary_", cv)), values_from = Count, values_fill = list(Count = 0)) %>%
    column_to_rownames("Macro_State") %>% as.matrix()
  if(nrow(tbl_data) >= 2 && ncol(tbl_data) >= 2 && sum(colSums(tbl_data) > 0) >= 2) {
    res <- fisher.test(tbl_data, simulate.p.value=TRUE)
    macro_stats <- macro_stats %>% add_row(Structural_Feature=cv, Target_State=TARGET_CATEGORIES[[cv]], p_value=res$p.value)
  }
}
if (nrow(macro_stats) > 0) write_csv(macro_stats %>% mutate(p.adj=p.adjust(p_value,"BH")) %>% add_significance("p.adj"), file.path(OUT_DIR, "Phase7_MacroState_Signatures_Stats.csv"))

master_out <- file.path(OUT_DIR, "Phase7_Complete_Structural_Metadata.csv")
write_csv(pca_metadata, master_out)
cat(sprintf("\n[✓] Unified Pipeline Complete! Master data saved to: %s\n", master_out))
