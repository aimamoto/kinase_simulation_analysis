# ==============================================================================
# SCRIPT: erbb_asymmetric_engine.R
# PURPOSE: Master structural, statistical, and ML evaluation for ERBB/Asymmetric 
#          Kinase AlphaFold3 ensembles.
#
# LOGICAL PIPELINE FLOW:
#   PART 1: Data Aggregation & Feature Engineering
#   PART 2: Global Thermodynamic Assembly (ATP-Agnostic Receiver Probability)
#   PART 3: 2ATP C-Spine Integrity Categorical Matrices
#   PART 4: 1ATP Competitive Binding & 2ATP Catalytic Clamping
#   PART 5: FAMILY-SPECIFIC Network Correlation & Density (MAC) (2ATP)
#   PART 6: FAMILY-SPECIFIC Unsupervised 3D Discovery (PCA & GMM) (2ATP)
#   PART 7: FAMILY-SPECIFIC KinCore Biological Signatures & Variant Composition
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
TARGET_TYPE <- if (length(args) > 0) args[1] else stop("Target Type required (e.g., EGFR).")
CLUSTER_METHOD <- if (length(args) > 1) tolower(args[2]) else "gmm"
FILE_LIST_PATH <- if (length(args) > 2) args[3] else stop("CSV list file required.")

if (!(CLUSTER_METHOD %in% c("kmeans", "gmm"))) {
  stop("[!] Invalid CLUSTER_METHOD. Must be 'kmeans' or 'gmm'.")
}

OUT_DIR <- sprintf("plots_and_stats_%s_%s_ERBB", TARGET_TYPE, toupper(CLUSTER_METHOD))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TARGET_CATEGORIES <- list(
  State = "Active (BLAminus)", C_Helix = "In", R_Spine = "Intact", Spatial = "DFGin"
)

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(rstatix)
  library(RVAideMemoire)
  library(patchwork)
  library(corrplot)
  library(factoextra)
  library(RColorBrewer)
  library(cluster)
  library(mclust)
  library(plotly)
  library(ggrepel)
})

cat(sprintf("[*] Asymmetric ERBB Engine Initialized.\n    Target: %s\n    Method: %s\n    Output: ./%s/\n", 
            TARGET_TYPE, toupper(CLUSTER_METHOD), OUT_DIR))

# --- THEMES & COLORS ---
PLOT_THEME <- theme_bw(base_size = 14) + 
  theme(text = element_text(color = "black"), axis.text = element_text(color = "black"),
        legend.position = "bottom", plot.title = element_text(face = "bold", size = 14))

ROLE_COLORS <- c("Activator" = "#F8766D", "Receiver" = "#00BFC4")

get_distinct_colors <- function(k) {
  tab10_bold <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  tab10_light <- c("#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
  tab20_smart <- c(tab10_bold, tab10_light)
  if (k <= 20) { return(tab20_smart[1:k]) } else { return(colorRampPalette(tab20_smart)(k)) }
}

# ==============================================================================
# PART 1: DATA AGGREGATION & FEATURE ENGINEERING
# ==============================================================================
cat("\n[*] Phase 1: Loading Data Provided by Python Wrapper...\n")

if(!file.exists(FILE_LIST_PATH)) stop(sprintf("[!] Cannot find file list at %s", FILE_LIST_PATH))
all_csvs <- readLines(FILE_LIST_PATH)
all_csvs <- all_csvs[all_csvs != ""]
if(length(all_csvs) == 0) stop("[!] No CSV files provided in the list.")

df_master <- map_dfr(all_csvs, ~read_csv(.x, col_types = cols(.default = "c"), na = c("", "NA", "N/A", "None"), show_col_types = FALSE)) %>%
  distinct(Simulation_ID, Chain, .keep_all = TRUE) %>%
  mutate(
    Receptor_A = str_match(Simulation_ID, "a-([^_]+)")[, 2],
    Receptor_B = str_match(Simulation_ID, "b-([^_]+)")[, 2],
    ATP_Count = str_match(Simulation_ID, "_(\\d+)atp")[, 2]
  ) %>%
  mutate(
    Receptor_A = str_remove_all(replace_na(Receptor_A, "Unknown"), "cattail"),
    Receptor_B = str_remove_all(replace_na(Receptor_B, "Unknown"), "cattail"),
    ATP_Count = replace_na(ATP_Count, "0"),
    Is_Homodimer = (Receptor_A == Receptor_B),
    Variant_A = ifelse(Receptor_A == tolower(TARGET_TYPE), "WT", toupper(str_remove(Receptor_A, paste0(tolower(TARGET_TYPE), "-")))),
    Variant_B = ifelse(Receptor_B == tolower(TARGET_TYPE), "WT", toupper(str_remove(Receptor_B, paste0(tolower(TARGET_TYPE), "-")))),
    Complex_Type = ifelse(Is_Homodimer, paste0(TARGET_TYPE, "_", Variant_A, "_HOMO"), paste(Receptor_A, "vs", Receptor_B, sep="_")),
    Plot_Group = case_when(Is_Homodimer ~ "All_Homodimers", TRUE ~ paste0(toupper(Receptor_A), "_Heterodimers"))
  )

df_master <- df_master %>%
  group_by(Simulation_ID) %>%
  arrange(Chain) %>%
  mutate(
    Chain_Index = row_number(),
    Chain_Identity = ifelse(Chain_Index == 1, Receptor_A, Receptor_B)
  ) %>%
  ungroup() %>%
  mutate(
    Clean_Identity = case_when(
      Chain_Identity == tolower(TARGET_TYPE) ~ "WT",
      str_detect(Chain_Identity, paste0(tolower(TARGET_TYPE), "-")) ~ toupper(str_remove(Chain_Identity, paste0(tolower(TARGET_TYPE), "-"))),
      TRUE ~ toupper(Chain_Identity)
    ),
    Is_Plot_Target = case_when(
      Is_Homodimer ~ TRUE,
      Receptor_A == tolower(TARGET_TYPE) & Receptor_B != tolower(TARGET_TYPE) & str_detect(Receptor_B, tolower(TARGET_TYPE)) ~ (Chain_Index == 2),
      Receptor_B == tolower(TARGET_TYPE) & Receptor_A != tolower(TARGET_TYPE) & str_detect(Receptor_A, tolower(TARGET_TYPE)) ~ (Chain_Index == 1),
      !str_detect(Receptor_A, tolower(TARGET_TYPE)) & str_detect(Receptor_B, tolower(TARGET_TYPE)) ~ (Chain_Index == 2),
      str_detect(Receptor_A, tolower(TARGET_TYPE)) & !str_detect(Receptor_B, tolower(TARGET_TYPE)) ~ (Chain_Index == 1),
      TRUE ~ (Chain_Index == 2)
    )
  )

dist_cols <- grep("_Dist$", colnames(df_master), value = TRUE)
df_master <- df_master %>% mutate(across(all_of(c("Phi_D", "Psi_D", dist_cols)), as.numeric))
df_master <- df_master %>% filter(Role %in% c("Activator", "Receiver"))
df_master$Role <- factor(df_master$Role, levels = c("Activator", "Receiver"))

# FIX: Explicitly unnest the list-columns generated by pivot_wider
sim_orientations <- df_master %>%
  select(Simulation_ID, Role, Clean_Identity) %>%
  group_by(Simulation_ID, Role) %>% slice(1) %>% ungroup() %>%
  pivot_wider(names_from = Role, values_from = Clean_Identity, values_fn = list) %>%
  mutate(
    Activator = sapply(Activator, function(x) if(is.null(x)) "Unk" else as.character(x[1])),
    Receiver = sapply(Receiver, function(x) if(is.null(x)) "Unk" else as.character(x[1]))
  ) %>%
  unnest(cols = c(Activator, Receiver)) %>%
  mutate(Dimer_Orientation = paste0(Activator, " (Act) ->\n", Receiver, " (Rec)")) %>% 
  select(Simulation_ID, Dimer_Orientation)

df_master <- df_master %>% left_join(sim_orientations, by = "Simulation_ID") %>%
  mutate(Dimer_Orientation = ifelse(Is_Homodimer, paste0(Clean_Identity, " (Act) ->\n", Clean_Identity, " (Rec)"), Dimer_Orientation),
         Group_Name = paste(Plot_Group, Complex_Type, Role, sep = "@"))

order_complexes_wt_first <- function(complex_names) {
  wt_homo <- complex_names[complex_names == paste0(TARGET_TYPE, "_WT_HOMO")]
  wt_hetero <- complex_names[str_detect(complex_names, paste0("vs_", tolower(TARGET_TYPE), "$")) & complex_names != paste0(TARGET_TYPE, "_WT_HOMO")]
  mut_names <- sort(complex_names[!complex_names %in% c(wt_homo, wt_hetero)])
  return(c(wt_homo, wt_hetero, mut_names))
}

cat(sprintf("    -> Parsed %d unique kinase chains.\n", nrow(df_master)))
cat("    -> NOTE: Asymmetric Engine detected. Bypassing Ligand Truth Overwrite.\n")

# ==============================================================================
# Helper Functions for Universal Statistical Outputs
# ==============================================================================
extract_pairwise <- function(df_pvals, has_role = TRUE) {
  if(has_role) {
    res <- df_pvals %>% separate(Group_1, into = c("Grp1", "Comp1", "Role1"), sep = "@", remove = FALSE, fill = "right") %>%
      separate(Group_2, into = c("Grp2", "Comp2", "Role2"), sep = "@", remove = FALSE, fill = "right")
  } else {
    res <- df_pvals %>% separate(Group_1, into = c("Grp1", "Comp1"), sep = "@", remove = FALSE, fill = "right") %>%
      separate(Group_2, into = c("Grp2", "Comp2"), sep = "@", remove = FALSE, fill = "right") %>%
      mutate(Role1 = "N/A", Role2 = "N/A")
  }
  res %>% mutate(
    Partner1 = case_when(Grp1 == "All_Homodimers" ~ TARGET_TYPE, TRUE ~ str_remove(Grp1, "_Heterodimers")),
    Partner2 = case_when(Grp2 == "All_Homodimers" ~ TARGET_TYPE, TRUE ~ str_remove(Grp2, "_Heterodimers")),
    Variant1 = case_when(str_detect(Comp1, "HOMO") | str_detect(Comp1, paste0("vs_", tolower(TARGET_TYPE), "$")) ~ "WT", TRUE ~ toupper(str_remove(Comp1, paste0(".*vs_", tolower(TARGET_TYPE), "-")))),
    Variant2 = case_when(str_detect(Comp2, "HOMO") | str_detect(Comp2, paste0("vs_", tolower(TARGET_TYPE), "$")) ~ "WT", TRUE ~ toupper(str_remove(Comp2, paste0(".*vs_", tolower(TARGET_TYPE), "-"))))
  )
}

format_fisher_pw <- function(fisher_raw, has_role = TRUE) {
  as.data.frame(as.table(fisher_raw$p.value)) %>% filter(!is.na(Freq)) %>% rename(Group_1=Var1, Group_2=Var2, p_adj=Freq) %>%
    mutate(Significance = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", TRUE ~ "ns")) %>% extract_pairwise(has_role = has_role)
}

safe_write_csv <- function(df, path) { if (nrow(df) > 0) write_csv(df, path) }

write_structured_comparisons <- function(pw_df, out_dir, prefix, has_role = TRUE) {
  safe_write_csv(if(has_role) filter(pw_df, Partner1 == Partner2 & Variant1 != Variant2 & Role1 == Role2) else filter(pw_df, Partner1 == Partner2 & Variant1 != Variant2), file.path(out_dir, paste0(prefix, "_IntraDimerType.csv")))
  safe_write_csv(if(has_role) filter(pw_df, Variant1 == Variant2 & Partner1 != Partner2 & Role1 == Role2) else filter(pw_df, Variant1 == Variant2 & Partner1 != Partner2), file.path(out_dir, paste0(prefix, "_InterDimerType.csv")))
  if(has_role) safe_write_csv(filter(pw_df, Partner1 == Partner2 & Variant1 == Variant2 & Role1 != Role2), file.path(out_dir, paste0(prefix, "_IntraComplex_Role.csv")))
}

run_global_fisher <- function(matrix_data) {
  tryCatch({ res <- fisher.test(matrix_data, simulate.p.value = TRUE, B = 5000); return(res$p.value) }, error = function(e) return(NA))
}

run_safe_pairwise_fisher <- function(matrix_data) {
  tryCatch({
    mat_clean <- matrix_data[, colSums(matrix_data) > 0, drop = FALSE]
    if(ncol(mat_clean) < 2 || nrow(mat_clean) < 2) return(NULL)
    return(fisher.multcomp(mat_clean, p.method = "fdr"))
  }, error = function(e) return(NULL))
}

# ==============================================================================
# PART 2: GLOBAL THERMODYNAMIC RECEIVER PROBABILITY (ATP-AGNOSTIC)
# ==============================================================================
cat("\n[*] Phase 2: Global Receiver Probabilities (H1 & H2)...\n")

df_p2 <- df_master %>%
  filter(Is_Plot_Target == TRUE) %>%
  mutate(
    Dimer_Family = case_when(Is_Homodimer ~ "WT\nBASE", TRUE ~ paste0(toupper(str_extract(ifelse(Chain_Index==2, Receptor_A, Receptor_B), "^[a-zA-Z0-9]+")), "\nHETERO")),
    Display_Name = Clean_Identity
  ) %>%
  group_by(Plot_Group, Complex_Type, Dimer_Family, Display_Name) %>%
  summarise(True_Total = n(), Receiver = sum(Role == "Receiver", na.rm = TRUE), Activator = sum(Role == "Activator", na.rm = TRUE), Unclassified = True_Total - (Receiver + Activator), Receiver_Pct = (Receiver / True_Total) * 100, .groups = "drop") %>% filter(True_Total > 0)

if (nrow(df_p2) > 0) {
  df_p2$Dimer_Family <- factor(df_p2$Dimer_Family, levels = c("WT\nBASE", sort(unique(df_p2$Dimer_Family[df_p2$Dimer_Family != "WT\nBASE"]))))
  unique_vars <- unique(df_p2$Display_Name)
  df_p2$Display_Name <- factor(df_p2$Display_Name, levels = rev(c("WT", sort(unique_vars[unique_vars != "WT"]))))
  
  plot_p2 <- ggplot(df_p2, aes(x = Receiver_Pct, y = Display_Name)) +
    geom_col(fill = "#2C7FB8", width = 0.7, color = "black", alpha = 0.85) +
    geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.8, alpha = 0.7) +
    geom_text(aes(x = 1, label = ifelse(Unclassified > 0, paste0(Receiver, "/", True_Total, " (", Unclassified, " unk)"), paste0(Receiver, "/", True_Total)), color = as.character(Receiver_Pct < 20)), hjust = 0, fontface = "bold", size = 4) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    facet_grid(Dimer_Family ~ ., scales = "free_y", space = "free_y") + PLOT_THEME + scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
    labs(title = sprintf("Receiver Probability of the Target %s Subunit", TARGET_TYPE), x = "Probability of Acting as Receiver (%)", y = "Variant") +
    theme(strip.background = element_rect(fill = "grey90"), strip.text.y = element_text(face = "bold", angle = 0))
  suppressWarnings(ggsave(file.path(OUT_DIR, "Plot_P2_Global_Receiver_Probability.pdf"), plot_p2, width = 12, height = 14, dpi = 300))
  
  df_global_roles <- df_p2 %>% mutate(Group_Name = paste(Plot_Group, Complex_Type, sep="@"))
  p2_global_res <- df_global_roles %>% group_by(Plot_Group) %>% group_modify(~ {
    mat <- .x %>% select(Receiver, Non_Receiver=Activator) %>% as.matrix()
    if(nrow(mat) > 1) tibble(Global_p_value = run_global_fisher(mat)) else tibble(Global_p_value = NA)
  }) %>% mutate(Significance = case_when(Global_p_value < 0.001 ~ "***", Global_p_value < 0.05 ~ "*", TRUE ~ "ns"))
  safe_write_csv(p2_global_res, file.path(OUT_DIR, "Stats_P2_GlobalRole_GLOBAL_Fisher.csv"))
  
  role_mat <- df_global_roles %>% select(Group_Name, Receiver, Activator) %>% column_to_rownames("Group_Name") %>% as.matrix()
  role_pw_raw <- run_safe_pairwise_fisher(role_mat)
  if(!is.null(role_pw_raw)) write_structured_comparisons(format_fisher_pw(role_pw_raw, has_role = FALSE), OUT_DIR, "Stats_P2_GlobalRole_Fisher", has_role = FALSE)
}


# ==============================================================================
# PART 3: C-SPINE INTEGRITY (CATEGORICAL)
# ==============================================================================
cat("[*] Phase 3: C-Spine Integrity Matrices...\n")

plot_cspine <- function(df, target_role) {
  df_agg <- df %>% filter(State != "N/A" & State != "Unknown", Role == target_role) %>%
    mutate(Dimer_Family = toupper(str_extract(Receptor_A, "^[a-zA-Z0-9]+")), ATP_State = paste0(ATP_Count, " ATP")) %>%
    group_by(Dimer_Family, Complex_Type, ATP_State, C_Spine) %>% summarise(Count = n(), .groups = "drop") %>%
    group_by(Dimer_Family, Complex_Type, ATP_State) %>% mutate(ATP_Total = sum(Count), Pct = Count / ATP_Total) %>%
    group_by(Dimer_Family, Complex_Type) %>% mutate(Overall_Total = sum(Count), Label_with_N = paste0(Complex_Type, " (N=", Overall_Total, ")")) %>% ungroup() %>% filter(ATP_Total > 0) %>%
    mutate(C_Spine = factor(C_Spine, levels = c("Intact", "Ligand Distant", "No Ligand")), Dimer_Family = factor(Dimer_Family, levels = c(TARGET_TYPE, unique(Dimer_Family[Dimer_Family != TARGET_TYPE]))))
  
  df_agg$Complex_Type <- factor(df_agg$Complex_Type, levels = rev(order_complexes_wt_first(unique(df_agg$Complex_Type))))
  
  ggplot(df_agg, aes(x = Pct, y = Complex_Type, fill = C_Spine)) + geom_col(position = "fill", color = "black", linewidth = 0.3, alpha = 0.9) +
    geom_text(aes(label = ifelse(Pct > 0.05, Count, "")), position = position_fill(vjust = 0.5), color = "gray10", fontface = "bold", size = 3.5) +
    facet_grid(Dimer_Family ~ ATP_State, scales = "free_y", space = "free_y") + PLOT_THEME + scale_x_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = c("Intact" = "#8DA0CB", "Ligand Distant" = "#66C2A5", "No Ligand" = "#FC8D62")) +
    labs(title = paste("Catalytic Spine Integrity vs. ATP Availability -", target_role), x = "Proportion of Structural Ensemble", y = "Simulation Dimer Pair", fill = "C-Spine Status") +
    theme(strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold"))
}
suppressWarnings(ggsave(file.path(OUT_DIR, "Plot_P3_CSpine_Role_Receiver.pdf"), plot_cspine(df_master, "Receiver"), width = 14, height = 12, dpi = 300))
suppressWarnings(ggsave(file.path(OUT_DIR, "Plot_P3_CSpine_Role_Activator.pdf"), plot_cspine(df_master, "Activator"), width = 14, height = 12, dpi = 300))


# ==============================================================================
# PART 4, 5, 6, 7: FAMILY-SPECIFIC DEEP DIVES (2ATP)
# ==============================================================================
cat("\n[*] Phases 4-7: Family-Specific Analytics (Binding, Allostery, PCA, GMM)...\n")

plot_groups <- unique(df_master$Plot_Group)
df_universal_anchor <- df_master %>% filter(Complex_Type == paste0(TARGET_TYPE, "_WT_HOMO"))

for (grp in plot_groups) {
  cat(sprintf("    -> Processing Family: %s\n", grp))
  family_out_dir <- file.path(OUT_DIR, paste0("Family_", grp, "_2ATP"))
  dir.create(family_out_dir, showWarnings = FALSE)
  
  df_grp <- df_master %>% filter(Plot_Group == grp)
  if (grp != "All_Homodimers" && nrow(df_universal_anchor) > 0) df_grp <- bind_rows(df_universal_anchor, df_grp)
  df_grp$Complex_Type <- factor(df_grp$Complex_Type, levels = order_complexes_wt_first(unique(df_grp$Complex_Type)))
  
  # --- 4A: 1ATP Preference ---
  df_grp_1atp <- df_grp %>% filter(ATP_Count == "1")
  if (nrow(df_grp_1atp) > 0) {
    total_sims_grp <- df_grp_1atp %>% group_by(Complex_Type, Dimer_Orientation) %>% summarise(Total_Sims = n_distinct(Simulation_ID), .groups = "drop")
    expected_grid_grp <- df_grp_1atp %>% distinct(Complex_Type, Dimer_Orientation) %>% crossing(Role = factor(c("Activator", "Receiver"), levels = c("Activator", "Receiver")))
    df_grp_1atp_winners <- df_grp_1atp %>% filter(C_Spine != "No Ligand") %>% group_by(Simulation_ID) %>% slice_min(order_by = HRD_ATP_Dist, n = 1, with_ties = FALSE) %>% ungroup()
    
    if (nrow(df_grp_1atp_winners) > 0) {
      df_grp_1atp_summary <- df_grp_1atp_winners %>% group_by(Complex_Type, Dimer_Orientation, Role) %>% summarise(Bound_ATP = n(), .groups = "drop") %>% 
        right_join(expected_grid_grp, by = c("Complex_Type", "Dimer_Orientation", "Role")) %>% left_join(total_sims_grp, by = c("Complex_Type", "Dimer_Orientation")) %>% 
        mutate(Bound_ATP = replace_na(Bound_ATP, 0), Binding_Percentage = (Bound_ATP / Total_Sims) * 100)
      
      df_grp_1atp_summary$Dimer_Orientation <- factor(df_grp_1atp_summary$Dimer_Orientation)
      
      p1atp <- ggplot(df_grp_1atp_summary, aes(x = Dimer_Orientation, y = Binding_Percentage, fill = Role)) + 
        geom_col(position = position_dodge(width = 0.8, preserve = "single"), color = "black", alpha = 0.8) + 
        geom_text(aes(y = ifelse(Binding_Percentage >= 20, 3, Binding_Percentage + 3), label = ifelse(Bound_ATP > 0, sprintf("%d\n(%.0f%%)", Bound_ATP, Binding_Percentage), "")), position = position_dodge(width = 0.8, preserve = "single"), vjust = 0, size = 3.5, fontface = "bold", color = "black") + 
        facet_wrap(~ Complex_Type, ncol = 3, scales = "free_x") + PLOT_THEME + scale_fill_manual(values = ROLE_COLORS, drop = FALSE) + 
        labs(title = sprintf("ATP Binding Preference (1ATP): %s", str_replace_all(grp, "_", " ")), x = "Dimer Configuration", y = "ATP Binding Success (%)") + 
        scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"), strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold"))
      suppressWarnings(ggsave(file.path(family_out_dir, sprintf("Plot_P4a_1ATP_Preference_%s.pdf", grp)), p1atp, width = 10, height = 10))
    }
  }
  
  # --- 4B/4C: 2ATP Local Statistics & Clamping Plots ---
  df_2atp_local <- df_grp %>% filter(ATP_Count == "2" & !is.na(HRD_ATP_Dist)) %>% group_by(Group_Name) %>% filter(n() >= 3) %>% ungroup()
  if (nrow(df_2atp_local) > 0 && n_distinct(df_2atp_local$Group_Name) >= 2) {
    tryCatch({
      if (sum(!is.na(df_2atp_local$HRD_ATP_Dist)) > 0) {
        p4b_pw_raw <- df_2atp_local %>% filter(!is.na(HRD_ATP_Dist)) %>% wilcox_test(HRD_ATP_Dist ~ Group_Name, p.adjust.method = "fdr") %>% add_significance("p.adj")
        write_structured_comparisons(extract_pairwise(select(p4b_pw_raw, Group_1=group1, Group_2=group2, p_adj=p.adj, Significance=p.adj.signif)), family_out_dir, "Stats_P4b_2ATP_Clamp_Wilcox")
      }
      if (sum(!is.na(df_2atp_local$aCb4_aE_Dist)) > 0) {
        p4c_pw_raw <- df_2atp_local %>% filter(!is.na(aCb4_aE_Dist)) %>% wilcox_test(aCb4_aE_Dist ~ Group_Name, p.adjust.method = "fdr") %>% add_significance("p.adj")
        write_structured_comparisons(extract_pairwise(select(p4c_pw_raw, Group_1=group1, Group_2=group2, p_adj=p.adj, Significance=p.adj.signif)), family_out_dir, "Stats_P4c_Allo_aCb4_Wilcox")
      }
    }, error = function(e) {})
  }
  
  plot_allo_metric <- function(df, metric_col, y_label, title_suffix) {
    df_filt <- df %>% filter(ATP_Count == "2" & !is.na(!!sym(metric_col)))
    if(nrow(df_filt) == 0) return(NULL)
    
    fill_args <- setNames(list(NA), metric_col)
    
    df_filt <- df_filt %>%
      complete(nesting(Complex_Type, Dimer_Orientation), Role, fill = fill_args) %>%
      filter(!is.na(Complex_Type))
    
    df_filt$Role <- factor(df_filt$Role, levels = c("Activator", "Receiver"))
    
    ggplot(df_filt, aes(x = Dimer_Orientation, y = !!sym(metric_col), fill = Role)) + 
      geom_violin(trim = FALSE, alpha = 0.6, color = NA, position = position_dodge(width = 0.8, preserve = "single")) + 
      geom_boxplot(aes(group = interaction(Dimer_Orientation, Role)), width = 0.25, fill = "white", color = "black", outlier.alpha = 0.3, position = position_dodge(width = 0.8, preserve = "single")) + 
      facet_wrap(~ Complex_Type, ncol = 3, scales = "free_x") + PLOT_THEME + scale_fill_manual(values = ROLE_COLORS, drop = FALSE) + 
      labs(title = paste("Structural Shift:", title_suffix, "-", str_replace_all(grp, "_", " ")), x = "Dimer Configuration", y = y_label) + 
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face="bold"), strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold"))
  }
  
  p2atp <- plot_allo_metric(df_grp, "HRD_ATP_Dist", "HRD-Asp to ATP Distance (Å)", "Catalytic Clamping")
  if(!is.null(p2atp)) suppressWarnings(ggsave(file.path(family_out_dir, sprintf("Plot_P4b_2ATP_Clamping_%s.pdf", grp)), p2atp + geom_hline(yintercept = 4.5, linetype="dashed", color="red"), width = 10, height = 10))
  
  pacb4 <- plot_allo_metric(df_grp, "aCb4_aE_Dist", "aC-b4 to aE Dist (Å)", "N/C-lobe Anchor")
  if(!is.null(pacb4)) suppressWarnings(ggsave(file.path(family_out_dir, sprintf("Plot_P4c_Allo_aCb4_%s.pdf", grp)), pacb4, width = 10, height = 10))
  
  
  # --- 5: NETWORK CORRELATION & MAC ---
  df_pca_base <- df_grp %>% filter(ATP_Count == "2")
  all_dist_cols <- grep("_Dist$", colnames(df_pca_base), value = TRUE)
  valid_dist_cols <- names(which(sapply(df_pca_base %>% select(all_of(all_dist_cols)), function(x) sum(!is.na(x)) >= 5)))
  
  if(length(valid_dist_cols) >= 2) {
    df_pca_complete <- df_pca_base %>% drop_na(all_of(valid_dist_cols))
    
    if(nrow(df_pca_complete) > 0) {
      high_contrast_pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)
      density_stats <- tibble(Group_Name = character(), Role = character(), N = numeric(), MAC = numeric())
      edge_weight_dists <- list() 
      
      pdf(file.path(family_out_dir, "Plot_P5_Correlation_Heatmaps_by_Role.pdf"), width = 28, height = 7 * ceiling(length(unique(df_pca_complete$Group_Name))/4))
      par(mfrow = c(ceiling(length(unique(df_pca_complete$Group_Name))/4), 4))
      
      for (sub_grp in unique(df_pca_complete$Group_Name)) {
        df_sub <- df_pca_complete %>% filter(Group_Name == sub_grp) %>% select(all_of(valid_dist_cols))
        if (nrow(df_sub) >= 5) {
          c_mat <- suppressWarnings(cor(df_sub, method = "spearman"))
          c_mat[is.na(c_mat)] <- 0 
          corrplot(c_mat, method="color", type="upper", tl.col="black", tl.cex=1.2, cl.cex=1.2, addgrid.col="white", title=sprintf("%s\n(n=%d)", str_replace(sub_grp, "@", "\n"), nrow(df_sub)), mar=c(0,0,4,0), col=high_contrast_pal)
          
          abs_edges <- abs(c_mat[upper.tri(c_mat)])
          edge_weight_dists[[sub_grp]] <- abs_edges
          density_stats <- density_stats %>% add_row(Group_Name=sub_grp, Role=str_split(sub_grp, "@")[[1]][3], N=nrow(df_sub), MAC=mean(abs_edges))
        }
      }
      dev.off()
      
      if (nrow(density_stats) > 0) {
        write_csv(density_stats %>% arrange(desc(MAC)), file.path(family_out_dir, "Stats_P5_Network_Density_MAC.csv"))
        p_mac <- density_stats %>%
          ggplot(aes(x = reorder(Group_Name, MAC), y = MAC, fill = Role)) + geom_col(color = "black", width = 0.7, alpha = 0.85) +
          coord_flip() + scale_fill_manual(values = ROLE_COLORS) + PLOT_THEME + labs(title = sprintf("%s: Global Network Rigidity (2ATP)", grp), subtitle = "Mean Absolute Correlation (MAC) of internal dynamics", x = NULL, y = "Global Coupling Score (MAC)")
        suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P5_Global_Network_Density.pdf"), p_mac, width = 12, height = 8))
        
        mac_stats <- tibble(Condition_1 = character(), Condition_2 = character(), MAC_1 = numeric(), MAC_2 = numeric(), p_value = numeric())
        if(length(edge_weight_dists) >= 2) {
          mac_pairs <- combn(names(edge_weight_dists), 2, simplify = FALSE)
          for(mp in mac_pairs) {
            p_val <- tryCatch({
              suppressWarnings(wilcox.test(edge_weight_dists[[mp[1]]], edge_weight_dists[[mp[2]]])$p.value)
            }, error = function(e) NA)
            
            mac_stats <- mac_stats %>% add_row(Condition_1 = mp[1], Condition_2 = mp[2], MAC_1 = mean(edge_weight_dists[[mp[1]]]), MAC_2 = mean(edge_weight_dists[[mp[2]]]), p_value = p_val)
          }
          
          mac_stats <- mac_stats %>% filter(!is.na(p_value))
          if(nrow(mac_stats) > 0) {
            mac_stats <- mac_stats %>% mutate(p.adj = p.adjust(p_value, "BH")) %>% add_significance("p.adj") %>% arrange(p.adj)
            write_csv(mac_stats, file.path(family_out_dir, "Stats_P5_Network_Density_Wilcox.csv"))
          }
        }
      }
      
      # --- 6: PCA & GMM CLUSTERING ---
      pca_data <- df_pca_complete %>% select(all_of(valid_dist_cols))
      pca_res <- prcomp(pca_data, scale. = TRUE) 
      
      if (CLUSTER_METHOD == "kmeans") {
        cluster_data <- pca_res$x[, 1:min(3, ncol(pca_res$x))]
        set.seed(42)
        gap_stat <- clusGap(cluster_data, FUN = kmeans, nstart = 25, K.max = 8, B = 100)
        optimal_k <- max(2, maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method="globalSEmax"))
        km_final <- kmeans(cluster_data, centers = optimal_k, nstart = 50)
        df_pca_complete$Macro_State <- factor(paste("State", km_final$cluster))
        
        suppressWarnings({
          p_gap <- fviz_gap_stat(gap_stat, maxSE=list(method="globalSEmax")) + PLOT_THEME
          ggsave(file.path(family_out_dir, "Plot_P6_Statistical_Proof_Gap.pdf"), p_gap, width = 8, height = 6)
        })
        
      } else {
        cluster_data <- pca_res$x[, 1:min(4, ncol(pca_res$x))]
        set.seed(42)
        gmm_res <- Mclust(cluster_data)
        optimal_k <- gmm_res$G
        df_pca_complete$Macro_State <- factor(paste("State", gmm_res$classification))
      }
      
      # Micro-Cluster Merging Algorithm
      state_counts <- df_pca_complete %>% group_by(Macro_State) %>% summarise(n=n(), .groups="drop")
      threshold <- max(10, nrow(df_pca_complete) * 0.05)
      small_states <- state_counts %>% filter(n < threshold) %>% pull(Macro_State)
      
      if(length(small_states) > 0 && length(unique(df_pca_complete$Macro_State)) > 1) {
        centroids <- df_pca_complete %>% group_by(Macro_State) %>% summarise(across(all_of(valid_dist_cols), mean)) %>% column_to_rownames("Macro_State")
        for (ss in small_states) {
          valid_targets <- setdiff(rownames(centroids), small_states)
          if(length(valid_targets) == 0) break
          dists <- as.matrix(dist(centroids))
          nearest <- names(which.min(dists[ss, valid_targets]))
          df_pca_complete <- df_pca_complete %>% mutate(Macro_State = as.character(Macro_State), Macro_State = ifelse(Macro_State == ss, nearest, Macro_State))
        }
        df_pca_complete$Macro_State <- factor(df_pca_complete$Macro_State)
        optimal_k <- length(unique(df_pca_complete$Macro_State))
        cat(sprintf("      -> Consolidated micro-clusters. Final Meta-States: %d\n", optimal_k))
      }
      
      state_colors <- get_distinct_colors(optimal_k)
      names(state_colors) <- levels(df_pca_complete$Macro_State)
      
      plot_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], State = df_pca_complete$Macro_State)
      hulls <- plot_df %>% group_by(State) %>% slice(chull(PC1, PC2)) %>% ungroup()
      
      p_pca <- ggplot(plot_df, aes(x = PC1, y = PC2, color = State, fill = State)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_polygon(data = hulls, alpha = 0.2, linewidth = 0) +
        scale_color_manual(values = state_colors) + scale_fill_manual(values = state_colors) +
        labs(title = sprintf("%s Meta-States", grp)) + PLOT_THEME
      suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P6_State_Clusters_PCA.pdf"), p_pca, width = 10, height = 7))
      
      # RESTORE 3D PLOTLY
      AXIS_STYLE <- list(titlefont = list(size = 18, color = "black"), tickfont = list(size = 14, color = "black"))
      p_3d <- plot_ly(data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2], PC3=pca_res$x[,3], State=df_pca_complete$Macro_State), 
                      x=~PC1, y=~PC2, z=~PC3, color=~State, colors=state_colors, type='scatter3d', mode='markers', 
                      marker=list(size=4, opacity=0.8, line=list(width=0))) %>%
        plotly::layout(title = list(text = sprintf("<b>%s Interactive 3D Phase Space</b>", grp), font = list(size = 24, color = "black")),
                       scene = list(xaxis = c(list(title = 'PC1'), AXIS_STYLE), yaxis = c(list(title = 'PC2'), AXIS_STYLE), zaxis = c(list(title = 'PC3'), AXIS_STYLE)))
      suppressWarnings(htmlwidgets::saveWidget(p_3d, file.path(family_out_dir, "Plot_P6_Interactive_3D_Space.html")))
      
      # State vs Role Composition
      expected_roles <- expand.grid(Macro_State = levels(df_pca_complete$Macro_State), Role = factor(c("Activator", "Receiver"), levels = c("Activator", "Receiver")))
      p_role_state <- df_pca_complete %>% 
        group_by(Macro_State, Role) %>% summarise(n = n(), .groups="drop") %>% 
        right_join(expected_roles, by = c("Macro_State", "Role")) %>% mutate(n = replace_na(n, 0)) %>%
        ggplot(aes(x = Macro_State, y = n, fill = Role)) + 
        geom_bar(stat="identity", position=position_dodge(preserve="single"), color="black") + 
        scale_fill_manual(values=ROLE_COLORS, drop=FALSE) + PLOT_THEME + labs(title="Do Activators and Receivers occupy different states?", y="Count")
      suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P6_State_Role_Composition.pdf"), p_role_state, width=8, height=6))
      
      p_variant_state <- df_pca_complete %>%
        group_by(Macro_State, Role, Complex_Type) %>% summarise(n = n(), .groups="drop") %>%
        ggplot(aes(x = Macro_State, y = n, fill = Complex_Type)) +
        geom_bar(stat="identity", position="stack", color="black") +
        facet_wrap(~ Role, scales = "free_y") + PLOT_THEME + 
        labs(title = "Variant Composition of Discovered States", y = "Count", fill = "Complex Type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P6_State_Variant_Composition.pdf"), p_variant_state, width=12, height=6))
      
      # --- 7: LOCAL KINCORE BIOLOGICAL SIGNATURES ---
      for(cv in names(TARGET_CATEGORIES)) {
        df_pca_complete[[paste0("Binary_", cv)]] <- factor(if_else(df_pca_complete[[cv]]==TARGET_CATEGORIES[[cv]], TARGET_CATEGORIES[[cv]], "Other"), levels=c(TARGET_CATEGORIES[[cv]], "Other"))
      }
      
      plot_feat <- function(v_name, o_var) {
        df_pca_complete %>% group_by(Macro_State, !!sym(v_name)) %>% summarise(n = n(), .groups="drop") %>%
          group_by(Macro_State) %>% mutate(Percent=n/sum(n)*100) %>%
          ggplot(aes(x=Macro_State, y=Percent, fill=!!sym(v_name))) + geom_bar(stat="identity", color="black", width=0.6) + 
          scale_fill_manual(values=setNames(c("#2171b5", "#cccccc"), c(TARGET_CATEGORIES[[o_var]], "Other"))) +
          labs(title=sprintf("%s Identity", o_var), x=NULL, y="Proportion (%)", fill="Feature") + PLOT_THEME
      }
      
      p_kincore <- (plot_feat("Binary_State","State") | plot_feat("Binary_C_Helix","C_Helix")) / 
        (plot_feat("Binary_R_Spine","R_Spine") | plot_feat("Binary_Spatial","Spatial"))
      suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P7_MacroState_Signatures.pdf"), p_kincore, width=16, height=12))
      
      macro_stats <- tibble(Structural_Feature=character(), Target_State=character(), p_value=numeric())
      for(cv in names(TARGET_CATEGORIES)) {
        tbl_data <- df_pca_complete %>% 
          group_by(Macro_State, !!sym(paste0("Binary_", cv))) %>% summarise(Count = n(), .groups="drop") %>%
          pivot_wider(names_from = !!sym(paste0("Binary_", cv)), values_from = Count, values_fill = list(Count = 0)) %>%
          column_to_rownames("Macro_State") %>% as.matrix()
        
        if(nrow(tbl_data) >= 2 && ncol(tbl_data) >= 2 && sum(colSums(tbl_data) > 0) >= 2) {
          res <- fisher.test(tbl_data, simulate.p.value=TRUE)
          macro_stats <- macro_stats %>% add_row(Structural_Feature=cv, Target_State=TARGET_CATEGORIES[[cv]], p_value=res$p.value)
        }
      }
      
      if (nrow(macro_stats) > 0) {
        write_csv(macro_stats %>% mutate(p.adj=p.adjust(p_value,"BH")) %>% add_significance("p.adj"), file.path(family_out_dir, "Stats_P7_MacroState_Signatures_Fisher.csv"))
      }
      
      write_csv(df_pca_complete, file.path(family_out_dir, "Phase7_Complete_Structural_Metadata.csv"))
    }
  } else {
    cat(sprintf("      [!] Not enough continuous data for PCA in family %s.\n", grp))
  }
}

cat(sprintf("\n[✓] ERBB Engine Complete! Nested outputs ready for Phase 8 Post-Hoc.\n"))