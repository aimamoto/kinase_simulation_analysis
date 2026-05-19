# ==============================================================================
# SCRIPT: erbb_asymmetric_engine.R
# PURPOSE: Master structural, statistical, and ML evaluation for ERBB/Asymmetric 
#          Kinase AlphaFold3 ensembles.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
TARGET_TYPE <- if (length(args) > 0) args[1] else stop("Target Type required.")
CLUSTER_METHOD <- if (length(args) > 1) tolower(args[2]) else "gmm"
FILE_LIST_PATH <- if (length(args) > 2) args[3] else stop("CSV list file required.")

if (!(CLUSTER_METHOD %in% c("kmeans", "gmm"))) {
  stop("[!] Invalid CLUSTER_METHOD. Must be 'kmeans' or 'gmm'.")
}

custom_out <- Sys.getenv("CUSTOM_OUT_DIR")
if (custom_out != "") { OUT_DIR <- custom_out } else {
  OUT_DIR <- sprintf("plots_and_stats_%s_%s_ERBB", TARGET_TYPE, toupper(CLUSTER_METHOD))
}
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

# ==============================================================================
# PLOT CONFIGURATION
# ==============================================================================
CFG <- list(
  font_size = 14,             
  title_size = 16,            
  canvas_w = 12,              
  canvas_h = 10,              
  p2_height = 14,             
  p3_h_scale = 0.6,           
  p4_h_scale = 2.5            
)

PLOT_THEME <- theme_bw(base_size = CFG$font_size) + 
  theme(text = element_text(color = "black"), axis.text = element_text(color = "black"),
        legend.position = "bottom", plot.title = element_text(face = "bold", size = CFG$title_size))

ROLE_COLORS <- c("Activator" = "#F8766D", "Receiver" = "#00BFC4")

get_distinct_colors <- function(k) {
  tab10_bold <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  tab10_light <- c("#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
  tab20_smart <- c(tab10_bold, tab10_light)
  if (k <= 20) return(tab20_smart[1:k]) else return(colorRampPalette(tab20_smart)(k))
}

# ==============================================================================
# PART 1: DATA AGGREGATION & LOCAL LIGAND STATE TRACKING
# ==============================================================================
cat("\n[*] Phase 1: Data Parsing & Explicit Identity Mapping...\n")

if(!file.exists(FILE_LIST_PATH)) stop(sprintf("[!] Cannot find file list at %s", FILE_LIST_PATH))
all_csvs <- readLines(FILE_LIST_PATH)
all_csvs <- all_csvs[all_csvs != ""]

df_master <- map_dfr(all_csvs, ~read_csv(.x, col_types = cols(.default = "c"), na = c("", "NA", "N/A", "None"), show_col_types = FALSE)) %>%
  distinct(Simulation_ID, Chain, .keep_all = TRUE)

df_types <- df_master %>% group_by(Simulation_ID) %>% arrange(Chain) %>%
  summarise(Raw_Type_A = Type[1], Raw_Type_B = if(n() > 1) Type[2] else Type[1], .groups = "drop") %>%
  mutate(Receptor_A = tolower(Raw_Type_A), Receptor_B = tolower(Raw_Type_B))

dist_cols <- grep("_Dist$", colnames(df_master), value = TRUE)
df_master <- df_master %>% mutate(across(all_of(c("Phi_D", "Psi_D", dist_cols)), as.numeric))

df_master <- df_master %>%
  left_join(df_types, by = "Simulation_ID") %>%
  mutate(
    ATP_Count = str_match(Simulation_ID, "_(\\d+)atp")[, 2],
    ATP_Count = as.numeric(replace_na(ATP_Count, "0")),
    Receptor_A = str_replace_all(Receptor_A, "(?i)cattail(?=-|$)|cat(?=-|$)", ""),
    Receptor_B = str_replace_all(Receptor_B, "(?i)cattail(?=-|$)|cat(?=-|$)", ""),
    Is_Homodimer = (Receptor_A == Receptor_B),
    Base_A = toupper(str_extract(Receptor_A, "^[a-zA-Z0-9]+")),
    Base_B = toupper(str_extract(Receptor_B, "^[a-zA-Z0-9]+")),
    Plot_Group = case_when(Is_Homodimer ~ "All_Homodimers", TRUE ~ paste0(Base_A, "_Heterodimers"))
  ) %>%
  group_by(Simulation_ID) %>% arrange(Chain) %>%
  mutate(
    Chain_Index = row_number(),
    Chain_Identity = ifelse(Chain_Index == 1, Receptor_A, Receptor_B),
    Base_Identity = ifelse(Chain_Index == 1, Base_A, Base_B),
    Clean_Identity = case_when(!str_detect(Chain_Identity, "-") ~ "WT", TRUE ~ toupper(str_remove(Chain_Identity, paste0(tolower(Base_Identity), "-")))),
    Complex_Type = case_when(Is_Homodimer ~ paste0(Base_Identity, "_", Clean_Identity, "_HOMO"), TRUE ~ paste(toupper(Receptor_A), "vs", toupper(Receptor_B), sep="_")),
    
    Ligand_State = case_when(
      ATP_Count >= 2 ~ "Holo",
      ATP_Count == 0 ~ "Apo",
      ATP_Count == 1 ~ {
        ranks <- rank(HRD_ATP_Dist, na.last = "keep", ties.method = "first")
        case_when(
          is.na(ranks) ~ "Apo",
          ranks %in% 1 & C_Spine != "No Ligand" ~ "Holo",
          TRUE ~ "Apo"
        )
      },
      TRUE ~ "Unknown"
    ),
    Cond_State = paste0(ATP_Count, "ATP_", Ligand_State)
  ) %>% ungroup()

df_master <- df_master %>% filter(Role %in% c("Activator", "Receiver")) %>% mutate(Role = factor(Role, levels = c("Activator", "Receiver")))

sim_orientations <- df_master %>% select(Simulation_ID, Role, Chain_Identity) %>%
  group_by(Simulation_ID, Role) %>% slice(1) %>% ungroup() %>%
  pivot_wider(names_from = Role, values_from = Chain_Identity, values_fn = list) %>%
  mutate(
    Activator = sapply(Activator, function(x) if(is.null(x)) "Unk" else toupper(as.character(x[1]))),
    Receiver = sapply(Receiver, function(x) if(is.null(x)) "Unk" else toupper(as.character(x[1])))
  ) %>%
  unnest(cols = c(Activator, Receiver)) %>% 
  mutate(
    Dimer_Orientation = paste0(Activator, " (Act) ->\n", Receiver, " (Rec)"),
    Orient_Clean = paste0(Activator, "->", Receiver)
  ) %>% select(Simulation_ID, Dimer_Orientation, Orient_Clean)

df_master <- df_master %>% left_join(sim_orientations, by = "Simulation_ID") %>%
  mutate(Group_Name = paste(Plot_Group, Complex_Type, Orient_Clean, Role, Cond_State, sep = "@"))

order_complexes_wt_first <- function(complex_names) {
  wt_homo <- complex_names[str_detect(complex_names, "_WT_HOMO$")]
  wt_hetero <- complex_names[str_detect(complex_names, "vs_[A-Z0-9]+$")]
  return(c(wt_homo, wt_hetero, sort(complex_names[!complex_names %in% c(wt_homo, wt_hetero)])))
}

# ==============================================================================
# HUMAN-READABLE STATISTICAL PARSING ENGINE (Upgraded Order & Polish)
# ==============================================================================
safe_write_csv <- function(df, path) { if (nrow(df) > 0) write_csv(df, path) }

format_pairwise_stats <- function(res_obj, is_directional = FALSE, metric_name = "Count") {
  is_rstatix <- "group1" %in% colnames(res_obj) || "group1" %in% names(res_obj)

  if (is_rstatix) {
    df <- as.data.frame(res_obj) %>% rename(Group_1 = group1, Group_2 = group2)
    if("p" %in% colnames(df)) df <- df %>% rename(p_raw = p) else df$p_raw <- NA
    if("p.adj" %in% colnames(df)) df <- df %>% rename(p_adj = p.adj) else {
       if("p_raw" %in% colnames(df) && !all(is.na(df$p_raw))) df$p_adj <- p.adjust(df$p_raw, method = "BH") else df$p_adj <- NA
    }
    if(".y." %in% colnames(df)) df <- df %>% rename(Metric = .y.) else df$Metric <- metric_name
    if("statistic" %in% colnames(df)) df <- df %>% rename(Test_Statistic = statistic) else df$Test_Statistic <- NA
    if("n1" %in% colnames(df)) df <- df %>% rename(N_1 = n1) else df$N_1 <- NA
    if("n2" %in% colnames(df)) df <- df %>% rename(N_2 = n2) else df$N_2 <- NA
  } else {
    # Matrix (Fisher) Handling
    df <- as.data.frame(as.table(res_obj$p.value)) %>% filter(!is.na(Freq)) %>% rename(Group_1 = Var1, Group_2 = Var2, p_adj = Freq) %>%
      mutate(Metric = metric_name, p_raw = NA, Test_Statistic = NA, N_1 = NA, N_2 = NA)
  }

  df <- df %>%
    mutate(Significance = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", TRUE ~ "ns")) %>%
    select(-any_of(c("p.adj.signif", "p.signif"))) 

  if (is_directional) {
    res <- df %>%
      separate(Group_1, into = c("Family1", "Dimer_Type_1", "Orientation_1", "Role_1", "Cond_State_1"), sep = "@", remove = FALSE, fill = "right", extra = "drop") %>%
      separate(Group_2, into = c("Family2", "Dimer_Type_2", "Orientation_2", "Role_2", "Cond_State_2"), sep = "@", remove = FALSE, fill = "right", extra = "drop") %>%
      select(any_of(c("Metric", "Dimer_Type_1", "Orientation_1", "Role_1", "Cond_State_1", "N_1",
                      "Dimer_Type_2", "Orientation_2", "Role_2", "Cond_State_2", "N_2",
                      "Test_Statistic", "p_raw", "p_adj", "Significance")),
             everything(), -Family1, -Family2, -Group_1, -Group_2)
  } else {
    res <- df %>%
      separate(Group_1, into = c("Family1", "Dimer_Type_1"), sep = "@", remove = FALSE, fill = "right", extra = "drop") %>%
      separate(Group_2, into = c("Family2", "Dimer_Type_2"), sep = "@", remove = FALSE, fill = "right", extra = "drop") %>%
      select(any_of(c("Metric", "Dimer_Type_1", "N_1", "Dimer_Type_2", "N_2",
                      "Test_Statistic", "p_raw", "p_adj", "Significance")),
             everything(), -Family1, -Family2, -Group_1, -Group_2)
  }

  # Clean up empty stat cols if purely non-parametric Fisher
  if(all(is.na(res$N_1))) res <- res %>% select(-N_1)
  if(all(is.na(res$N_2))) res <- res %>% select(-N_2)
  if(all(is.na(res$Test_Statistic))) res <- res %>% select(-Test_Statistic)
  if(all(is.na(res$p_raw))) res <- res %>% select(-p_raw)

  return(res)
}

write_stats_csvs <- function(pw_df, out_dir, prefix, is_directional = FALSE) {
  safe_write_csv(pw_df, file.path(out_dir, paste0(prefix, "_ALL_Comparisons.csv")))
  if (is_directional) {
    safe_write_csv(filter(pw_df, Dimer_Type_1 == Dimer_Type_2 & Orientation_1 == Orientation_2 & Role_1 != Role_2), file.path(out_dir, paste0(prefix, "_IntraDimer_RoleEffect.csv")))
    safe_write_csv(filter(pw_df, Dimer_Type_1 == Dimer_Type_2 & Orientation_1 != Orientation_2 & Role_1 == Role_2), file.path(out_dir, paste0(prefix, "_IntraDimer_DirectionalityEffect.csv")))
    safe_write_csv(filter(pw_df, Dimer_Type_1 != Dimer_Type_2 & Role_1 == Role_2), file.path(out_dir, paste0(prefix, "_InterDimer_Variants.csv")))
  }
}

run_safe_pairwise_fisher <- function(matrix_data) {
  tryCatch({
    mat_clean <- matrix_data[, colSums(matrix_data) > 0, drop = FALSE]
    if(ncol(mat_clean) < 2 || nrow(mat_clean) < 2) return(NULL)
    return(fisher.multcomp(mat_clean, p.method = "fdr"))
  }, error = function(e) return(NULL))
}

# ==============================================================================
# FAMILY-SPECIFIC LOOP (Phases 2 through 7)
# ==============================================================================
cat("\n[*] Executing Local Deep Dives (Phases 2-7)...\n")

plot_groups <- unique(df_master$Plot_Group)
df_universal_anchor <- df_master %>% filter(Plot_Group == "All_Homodimers")
ultimate_pca_master_out <- data.frame()

for (grp in plot_groups) {
  if(grp == "All_Homodimers" && length(plot_groups) > 1) next 
  cat(sprintf("\n    >>> Processing Family: %s\n", grp))
  
  family_out_dir <- file.path(OUT_DIR, paste0("Family_", grp))
  dir.create(family_out_dir, showWarnings = FALSE)
  
  base_identity <- str_remove(grp, "_Heterodimers")
  if(base_identity == "All_Homodimers") {
      df_grp <- df_master %>% filter(Plot_Group == grp)
  } else {
      relevant_homos <- df_universal_anchor %>% filter(Base_A == base_identity)
      df_grp <- bind_rows(relevant_homos, df_master %>% filter(Plot_Group == grp))
  }
  
  if(nrow(df_grp) == 0) next
  df_grp$Complex_Type <- factor(df_grp$Complex_Type, levels = order_complexes_wt_first(unique(df_grp$Complex_Type)))
  
  # --------------------------------------------------------------------------
  # PHASE 2: LOCAL THERMODYNAMIC PROBABILITY 
  # --------------------------------------------------------------------------
  df_p2 <- df_grp %>%
    filter(Base_Identity == base_identity) %>%
    mutate(
      Partner_Identity = ifelse(Chain_Index == 1, Receptor_B, Receptor_A),
      Partner_Base = toupper(str_extract(Partner_Identity, "^[a-zA-Z0-9]+")),
      Clean_Partner = case_when(!str_detect(Partner_Identity, "-") ~ "WT", TRUE ~ toupper(str_remove(Partner_Identity, paste0(tolower(Partner_Base), "-")))),
      Dimer_Family = case_when(Is_Homodimer ~ "WT\nBASE", TRUE ~ paste0(Partner_Base, "\nHETERO")),
      Display_Name = case_when(Is_Homodimer ~ paste(Base_Identity, Clean_Identity), TRUE ~ paste0(Base_Identity, " ", Clean_Identity, "  [+ ", Partner_Base, " ", Clean_Partner, "]"))
    ) %>%
    group_by(Plot_Group, Complex_Type, Dimer_Family, Display_Name) %>%
    summarise(
      True_Total = n(), Receiver = sum(Role == "Receiver", na.rm = TRUE), Activator = sum(Role == "Activator", na.rm = TRUE), 
      Unclassified = True_Total - (Receiver + Activator), Receiver_Pct = (Receiver / True_Total) * 100, .groups = "drop"
    ) %>% filter(True_Total > 0)

  if (nrow(df_p2) > 0) {
    df_p2$Dimer_Family <- factor(df_p2$Dimer_Family, levels = c("WT\nBASE", sort(unique(df_p2$Dimer_Family[df_p2$Dimer_Family != "WT\nBASE"]))))
    df_p2 <- df_p2 %>% arrange(Dimer_Family, Display_Name)
    df_p2$Display_Name <- factor(df_p2$Display_Name, levels = rev(unique(df_p2$Display_Name)))
    
    plot_p2 <- ggplot(df_p2, aes(x = Receiver_Pct, y = Display_Name)) +
      geom_col(fill = "#2C7FB8", width = 0.7, color = "black", alpha = 0.85) +
      geom_vline(xintercept = 50, linetype = "dashed", color = "red", linewidth = 0.8, alpha = 0.7) +
      geom_text(aes(x = 1, label = ifelse(Unclassified > 0, paste0(Receiver, "/", True_Total, " (", Unclassified, " unk)"), paste0(Receiver, "/", True_Total)), color = as.character(Receiver_Pct < 20)), hjust = 0, fontface = "bold", size = 4) +
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
      facet_grid(Dimer_Family ~ ., scales = "free_y", space = "free_y") + PLOT_THEME + scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
      labs(title = sprintf("Receiver Probability of the Target %s Subunit", base_identity), x = "Probability of Acting as Receiver (%)", y = "Variant Configuration") +
      theme(strip.background = element_rect(fill = "grey90"), strip.text.y = element_text(face = "bold", angle = 0))
    suppressWarnings(ggsave(file.path(family_out_dir, sprintf("Plot_P2_Receiver_Probability_%s.pdf", base_identity)), plot_p2, width = CFG$canvas_w, height = CFG$p2_height, dpi = 300))
    
    df_global_roles <- df_p2 %>% mutate(Test_Group = paste(Plot_Group, Complex_Type, sep="@"))
    p2_global_res <- df_global_roles %>% group_by(Plot_Group) %>% group_modify(~ {
      mat <- .x %>% select(Receiver, Non_Receiver=Activator) %>% as.matrix()
      if(nrow(mat) > 1) tibble(Global_p_value = tryCatch(fisher.test(mat, simulate.p.value=TRUE, B=5000)$p.value, error=function(e) NA)) else tibble(Global_p_value = NA)
    }) %>% mutate(Significance = case_when(Global_p_value < 0.001 ~ "***", Global_p_value < 0.05 ~ "*", TRUE ~ "ns"))
    safe_write_csv(p2_global_res, file.path(family_out_dir, "Stats_P2_GlobalRole_GLOBAL_Fisher.csv"))
    
    role_mat <- df_global_roles %>% distinct(Test_Group, .keep_all = TRUE) %>% select(Test_Group, Receiver, Activator) %>% column_to_rownames("Test_Group") %>% as.matrix()
    role_pw_raw <- run_safe_pairwise_fisher(role_mat)
    if(!is.null(role_pw_raw)) write_stats_csvs(format_pairwise_stats(role_pw_raw, is_directional=FALSE), family_out_dir, "Stats_P2_GlobalRole_Fisher", is_directional=FALSE)
  }

  # --------------------------------------------------------------------------
  # PHASE 3: LOCAL C-SPINE INTEGRITY 
  # --------------------------------------------------------------------------
  calc_p3_stats_and_plot <- function(df_in, target_role) {
    df_agg <- df_in %>% filter(State != "N/A" & State != "Unknown", Role == target_role) %>%
      mutate(
        C_Spine_Bin = ifelse(C_Spine == "Intact", "Intact", "Broken/Distant"),
        Test_Group = paste(Plot_Group, Complex_Type, sep="@")
      )
    
    if(nrow(df_agg) == 0) return(NULL)
    
    tbl <- table(df_agg$Test_Group, df_agg$C_Spine_Bin)
    if(nrow(tbl) > 1 && ncol(tbl) > 1) {
      pval <- tryCatch(fisher.test(tbl, simulate.p.value=TRUE, B=5000)$p.value, error=function(e) NA)
      write_csv(tibble(Role=target_role, Global_p_value=pval) %>% mutate(Significance = case_when(Global_p_value < 0.001 ~ "***", Global_p_value < 0.01 ~ "**", Global_p_value < 0.05 ~ "*", TRUE ~ "ns")), file.path(family_out_dir, paste0("Stats_P3_CSpine_GLOBAL_", target_role, ".csv")))
      
      pw <- run_safe_pairwise_fisher(as.matrix(tbl))
      if(!is.null(pw)) {
        write_stats_csvs(format_pairwise_stats(pw, is_directional=FALSE), family_out_dir, paste0("Stats_P3_CSpine_", target_role), is_directional=FALSE)
      }
    }

    df_plot <- df_agg %>%
      group_by(Complex_Type, Cond_State, C_Spine) %>% summarise(Count = n(), .groups = "drop") %>%
      group_by(Complex_Type, Cond_State) %>% mutate(Cond_Total = sum(Count), Pct = Count / Cond_Total) %>%
      group_by(Complex_Type) %>% mutate(Overall_Total = sum(Count), Label_with_N = paste0(Complex_Type, " (N=", Overall_Total, ")")) %>% ungroup() %>% filter(Cond_Total > 0) %>%
      mutate(C_Spine = factor(C_Spine, levels = c("Intact", "Ligand Distant", "No Ligand")))
      
    df_plot$Complex_Type <- factor(df_plot$Complex_Type, levels = rev(order_complexes_wt_first(unique(df_plot$Complex_Type))))
    plot_height <- max(6, n_distinct(df_plot$Complex_Type) * CFG$p3_h_scale)
    
    p <- ggplot(df_plot, aes(x = Pct, y = Complex_Type, fill = C_Spine)) + geom_col(position = "fill", color = "black", linewidth = 0.3, alpha = 0.9) +
      geom_text(aes(label = ifelse(Pct > 0.05, Count, "")), position = position_fill(vjust = 0.5), color = "gray10", fontface = "bold", size = 3.5) +
      facet_grid(~ Cond_State, scales = "free_x", space = "free_x") + PLOT_THEME + scale_x_continuous(labels = percent_format(accuracy = 1)) +
      scale_fill_manual(values = c("Intact" = "#8DA0CB", "Ligand Distant" = "#66C2A5", "No Ligand" = "#FC8D62")) +
      labs(title = paste("Catalytic Spine Integrity vs. ATP/Ligand State -", target_role), x = "Proportion of Structural Ensemble", y = "Simulation Dimer Pair", fill = "C-Spine Status") +
      theme(strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold"))
      
    return(list(plot=p, height=plot_height))
  }
  
  res_p3_rec <- calc_p3_stats_and_plot(df_grp, "Receiver")
  if(!is.null(res_p3_rec)) suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P3_CSpine_Role_Receiver.pdf"), res_p3_rec$plot, width = CFG$canvas_w + 3, height = res_p3_rec$height, dpi = 300))
  res_p3_act <- calc_p3_stats_and_plot(df_grp, "Activator")
  if(!is.null(res_p3_act)) suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P3_CSpine_Role_Activator.pdf"), res_p3_act$plot, width = CFG$canvas_w + 3, height = res_p3_act$height, dpi = 300))

  # --------------------------------------------------------------------------
  # PHASE 4a: 1ATP BINDING (TUG-OF-WAR)
  # --------------------------------------------------------------------------
  df_grp_1atp <- df_grp %>% filter(ATP_Count == 1)
  if (nrow(df_grp_1atp) > 0) {
    total_sims_grp <- df_grp_1atp %>% group_by(Complex_Type, Orient_Clean) %>% summarise(Total_Sims = n_distinct(Simulation_ID), .groups = "drop")
    expected_grid_grp <- df_grp_1atp %>% distinct(Complex_Type, Orient_Clean) %>% crossing(Role = factor(c("Activator", "Receiver"), levels = c("Activator", "Receiver")))
    
    df_grp_1atp_winners <- df_grp_1atp %>% filter(Ligand_State == "Holo") 
    
    if (nrow(df_grp_1atp_winners) > 0) {
      df_grp_1atp_summary <- df_grp_1atp_winners %>% group_by(Complex_Type, Orient_Clean, Role) %>% summarise(Bound_ATP = n(), .groups = "drop") %>% 
        right_join(expected_grid_grp, by = c("Complex_Type", "Orient_Clean", "Role")) %>% left_join(total_sims_grp, by = c("Complex_Type", "Orient_Clean")) %>% 
        mutate(Bound_ATP = replace_na(Bound_ATP, 0), Binding_Percentage = (Bound_ATP / Total_Sims) * 100)
      
      p4a_stats <- df_grp_1atp_summary %>%
        group_by(Complex_Type, Orient_Clean) %>%
        summarise(Total_Bound = sum(Bound_ATP), Act_Bound = sum(Bound_ATP[Role == "Activator"]), Rec_Bound = sum(Bound_ATP[Role == "Receiver"]), .groups = "drop") %>%
        filter(Total_Bound > 0) %>% rowwise() %>%
        mutate(p_value = binom.test(Act_Bound, Total_Bound, p = 0.5)$p.value, Preference = case_when(Act_Bound > Rec_Bound ~ "Activator", Act_Bound < Rec_Bound ~ "Receiver", TRUE ~ "None")) %>%
        ungroup() %>% mutate(p.adj = p.adjust(p_value, method = "BH")) %>% add_significance("p.adj")
      write_csv(p4a_stats, file.path(family_out_dir, "Stats_P4a_1ATP_Binding_Binomial.csv"))

      p4a_mat <- df_grp_1atp_summary %>%
        mutate(Test_Group = paste(grp, Complex_Type, sep="@")) %>%
        group_by(Test_Group) %>%
        summarise(Act_Bound = sum(Bound_ATP[Role == "Activator"]), Rec_Bound = sum(Bound_ATP[Role == "Receiver"]), .groups = "drop") %>%
        column_to_rownames("Test_Group") %>% as.matrix()

      p4a_pw <- run_safe_pairwise_fisher(p4a_mat)
      if(!is.null(p4a_pw)) write_stats_csvs(format_pairwise_stats(p4a_pw, is_directional=FALSE), family_out_dir, "Stats_P4a_1ATP_Binding_Fisher", is_directional=FALSE)

      plot_width <- max(CFG$canvas_w, n_distinct(df_grp_1atp_summary$Complex_Type) * 3)
      p1atp <- ggplot(df_grp_1atp_summary, aes(x = Orient_Clean, y = Binding_Percentage, fill = Role)) + 
        geom_col(position = position_dodge(width = 0.8, preserve = "single"), color = "black", alpha = 0.8) + 
        geom_text(aes(y = ifelse(Binding_Percentage >= 20, 3, Binding_Percentage + 3), label = ifelse(Bound_ATP > 0, sprintf("%d\n(%.0f%%)", Bound_ATP, Binding_Percentage), "")), position = position_dodge(width = 0.8, preserve = "single"), vjust = 0, size = 3.5, fontface = "bold", color = "black") + 
        facet_wrap(~ Complex_Type, ncol = 3, scales = "free_x") + PLOT_THEME + scale_fill_manual(values = ROLE_COLORS, drop = FALSE) + 
        labs(title = sprintf("ATP Binding Preference (1ATP): %s", str_replace_all(grp, "_", " ")), x = "Dimer Configuration", y = "ATP Binding Success (%)") + 
        scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"), strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold"))
      suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P4a_1ATP_Preference.pdf"), p1atp, width = plot_width/2, height = 12))
    }
  }
  
  # --------------------------------------------------------------------------
  # PHASE 4b/c: CONTINUOUS CLAMPING / ALLOSTERY
  # --------------------------------------------------------------------------
  df_local <- df_grp %>% filter(!is.na(HRD_ATP_Dist) | !is.na(aCb4_aE_Dist)) %>% group_by(Group_Name) %>% filter(n() >= 3) %>% ungroup()
  if (nrow(df_local) > 0 && n_distinct(df_local$Group_Name) >= 2) {
    tryCatch({
      if (sum(!is.na(df_local$HRD_ATP_Dist)) > 0) {
        p4b_pw_raw <- df_local %>% filter(!is.na(HRD_ATP_Dist)) %>% wilcox_test(HRD_ATP_Dist ~ Group_Name, p.adjust.method = "fdr")
        write_stats_csvs(format_pairwise_stats(p4b_pw_raw, is_directional=TRUE), family_out_dir, "Stats_P4b_Clamp_Wilcox", is_directional=TRUE)
      }
      if (sum(!is.na(df_local$aCb4_aE_Dist)) > 0) {
        p4c_pw_raw <- df_local %>% filter(!is.na(aCb4_aE_Dist)) %>% wilcox_test(aCb4_aE_Dist ~ Group_Name, p.adjust.method = "fdr")
        write_stats_csvs(format_pairwise_stats(p4c_pw_raw, is_directional=TRUE), family_out_dir, "Stats_P4c_Allo_aCb4_Wilcox", is_directional=TRUE)
      }
    }, error = function(e) {})
  }
  
  plot_allo_metric <- function(df_in, metric_col, y_label, title_suffix) {
    df_filt <- df_in %>% filter(!is.na(!!sym(metric_col)))
    if(nrow(df_filt) == 0) return(NULL)
    
    df_filt <- df_filt %>% complete(nesting(Complex_Type, Orient_Clean, Cond_State), Role, fill = setNames(list(NA), metric_col)) %>% filter(!is.na(Complex_Type))
    
    plot_width <- max(CFG$canvas_w, n_distinct(df_filt$Complex_Type) * 3)
    plot_height <- max(6, n_distinct(df_filt$Cond_State) * 3 + 3)
    
    p <- ggplot(df_filt, aes(x = Orient_Clean, y = !!sym(metric_col), fill = Role)) + 
      geom_violin(trim = FALSE, alpha = 0.6, color = NA, position = position_dodge(width = 0.8, preserve = "single")) + 
      geom_boxplot(aes(group = interaction(Orient_Clean, Role)), width = 0.25, fill = "white", color = "black", outlier.alpha = 0.3, position = position_dodge(width = 0.8, preserve = "single")) + 
      facet_grid(Cond_State ~ Complex_Type, scales = "free_x") + PLOT_THEME + scale_fill_manual(values = ROLE_COLORS, drop = FALSE) + 
      labs(title = paste("Structural Shift:", title_suffix, "-", str_replace_all(grp, "_", " ")), x = "Dimer Configuration", y = y_label) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face="bold"), strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold"))
    return(list(plot = p, width = plot_width, height = plot_height))
  }
  
  res_p4b <- plot_allo_metric(df_grp, "HRD_ATP_Dist", "HRD-Asp to ATP Distance (Å)", "Catalytic Clamping")
  if(!is.null(res_p4b)) suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P4b_Clamping.pdf"), res_p4b$plot + geom_hline(yintercept = 4.5, linetype="dashed", color="red"), width = res_p4b$width, height = res_p4b$height))
  res_p4c <- plot_allo_metric(df_grp, "aCb4_aE_Dist", "aC-b4 to aE Dist (Å)", "N/C-lobe Anchor")
  if(!is.null(res_p4c)) suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P4c_Allo_aCb4.pdf"), res_p4c$plot, width = res_p4c$width, height = res_p4c$height))
  
  # --------------------------------------------------------------------------
  # PHASE 5, 6, 7: UNSUPERVISED NETWORK ANALYSIS 
  # --------------------------------------------------------------------------
  df_pca_base <- df_grp 
  ligand_cols <- c("HRD_ATP_Dist", "DFG_Mg_Dist", "DFG_ATP_Dist", "PLoop_ATP_Dist", "Spine_Bridge_Dist")
  all_dist_cols <- setdiff(grep("_Dist$", colnames(df_pca_base), value = TRUE), ligand_cols)
  valid_dist_cols <- names(which(sapply(df_pca_base %>% select(all_of(all_dist_cols)), function(x) sum(!is.na(x)) >= 5)))
  
  if(length(valid_dist_cols) >= 2) {
    df_pca_complete <- df_pca_base %>% drop_na(all_of(valid_dist_cols))
    
    if(nrow(df_pca_complete) > 0) {
      
      df_pca_complete <- df_pca_complete %>% mutate(Plot_Network = paste0(Orient_Clean, " (", Role, " | ", Cond_State, ")"))
        
      high_contrast_pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)
      density_stats <- tibble(Stat_Group = character(), Plot_Network = character(), Role = character(), Cond_State = character(), N = numeric(), MAC = numeric())
      edge_weight_dists <- list() 
      
      pdf(file.path(family_out_dir, "Plot_P5_Correlation_Heatmaps.pdf"), width = 28, height = 7 * ceiling(length(unique(df_pca_complete$Group_Name))/4))
      par(mfrow = c(ceiling(length(unique(df_pca_complete$Group_Name))/4), 4))
      for (sub_grp in unique(df_pca_complete$Group_Name)) {
        df_sub <- df_pca_complete %>% filter(Group_Name == sub_grp) %>% select(all_of(valid_dist_cols))
        if (nrow(df_sub) >= 5) {
          c_mat <- suppressWarnings(cor(df_sub, method = "spearman"))
          c_mat[is.na(c_mat)] <- 0 
          plot_title_hm <- str_replace_all(unique(df_pca_complete$Plot_Network[df_pca_complete$Group_Name == sub_grp])[1], " \\| ", "\n")
          corrplot(c_mat, method="color", type="upper", tl.col="black", tl.cex=1.2, cl.cex=1.2, addgrid.col="white", title=sprintf("%s\n(n=%d)", plot_title_hm, nrow(df_sub)), mar=c(0,0,4,0), col=high_contrast_pal)
          abs_edges <- abs(c_mat[upper.tri(c_mat)])
          edge_weight_dists[[sub_grp]] <- abs_edges
          clean_plot_network <- unique(df_pca_complete$Plot_Network[df_pca_complete$Group_Name == sub_grp])[1]
          density_stats <- density_stats %>% add_row(Stat_Group=sub_grp, Plot_Network=clean_plot_network, Role=str_split(sub_grp, "@")[[1]][4], Cond_State=str_split(sub_grp, "@")[[1]][5], N=nrow(df_sub), MAC=mean(abs_edges))
        }
      }
      dev.off()
      
      if (nrow(density_stats) > 0) {
        write_csv(density_stats %>% select(-Plot_Network) %>% arrange(desc(MAC)), file.path(family_out_dir, "Stats_P5_Network_Density_MAC.csv"))
        
        p_mac <- density_stats %>% ggplot(aes(x = reorder(Plot_Network, MAC), y = MAC, fill = Role)) + geom_col(color = "black", width = 0.7, alpha = 0.85) + coord_flip() + scale_fill_manual(values = ROLE_COLORS) + PLOT_THEME + 
          labs(title = sprintf("%s: Global Network Rigidity", grp), subtitle = "Mean Absolute Correlation (MAC). Higher MAC = Greater Rigidity & Stronger Structural Coupling.", x = NULL, y = "Global Coupling Score (MAC)")
        suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P5_Global_Network_Density.pdf"), p_mac, width = CFG$canvas_w + 2, height = max(8, nrow(density_stats) * 0.4)))
        
        mac_stats <- tibble(group1 = character(), group2 = character(), MAC_1 = numeric(), MAC_2 = numeric(), p = numeric())
        if(length(edge_weight_dists) >= 2) {
          mac_pairs <- combn(names(edge_weight_dists), 2, simplify = FALSE)
          for(mp in mac_pairs) {
            p_val <- tryCatch(suppressWarnings(wilcox.test(edge_weight_dists[[mp[1]]], edge_weight_dists[[mp[2]]])$p.value), error = function(e) NA)
            mac_stats <- mac_stats %>% add_row(group1 = mp[1], group2 = mp[2], MAC_1 = mean(edge_weight_dists[[mp[1]]]), MAC_2 = mean(edge_weight_dists[[mp[2]]]), p = p_val)
          }
          mac_stats <- mac_stats %>% filter(!is.na(p)) %>% mutate(p.adj = p.adjust(p, "BH"), Delta_MAC = MAC_1 - MAC_2)
          if(nrow(mac_stats) > 0) {
            mac_pw <- format_pairwise_stats(mac_stats, is_directional=TRUE, metric_name="MAC")
            write_stats_csvs(mac_pw, family_out_dir, "Stats_P5_Network_Density_Wilcox", is_directional=TRUE)
            
            # --- VOLCANO PLOT BLOCK ---
            plot_df_volc <- mac_pw %>%
              mutate(
                LogP = -log10(p_adj + 1e-16),
                Significant = p_adj < 0.05,
                Comparison_Label = case_when(
                  p_adj < 0.01 & abs(Delta_MAC) > 0.03 ~ paste0(Dimer_Type_1, "(", substr(Role_1,1,3), ")\nvs\n", Dimer_Type_2, "(", substr(Role_2,1,3), ")"),
                  TRUE ~ ""
                )
              )
            
            p_volcano <- ggplot(plot_df_volc, aes(x = Delta_MAC, y = LogP)) +
              geom_point(aes(color = Significant), alpha = 0.8, size = 3) +
              geom_text_repel(aes(label = Comparison_Label), size = 3, max.overlaps = 15, box.padding = 0.5) +
              geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
              geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
              scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey80")) +
              PLOT_THEME + 
              labs(
                title = sprintf("Network Rigidity Shift (Volcano): %s", grp),
                subtitle = "Pairwise comparisons of MAC (Correlation) shifts",
                x = expression(Delta~"MAC (Group 1 - Group 2)"),
                y = expression("-log"[10]*"(FDR p-value)")
              ) + theme(legend.position = "none")
            
            suppressWarnings(ggsave(file.path(family_out_dir, "Plot_P5_Network_Density_Volcano.pdf"), p_volcano, width = CFG$canvas_w, height = 8))
          }
        }
      }
      
      df_pca_master_out <- data.frame()
      for (target_role in c("Activator", "Receiver")) {
        for (target_cond in unique(df_pca_complete$Cond_State)) {
          df_pca_role <- df_pca_complete %>% filter(Role == target_role, Cond_State == target_cond)
          if(nrow(df_pca_role) < 10) next
          
          pca_data <- df_pca_role %>% select(all_of(valid_dist_cols))
          pca_res <- prcomp(pca_data, scale. = TRUE) 
          
          if (CLUSTER_METHOD == "kmeans") {
            cluster_data <- pca_res$x[, 1:min(3, ncol(pca_res$x))]
            set.seed(42)
            gap_stat <- clusGap(cluster_data, FUN = kmeans, nstart = 25, K.max = 8, B = 100)
            optimal_k <- max(2, maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method="globalSEmax"))
            km_final <- kmeans(cluster_data, centers = optimal_k, nstart = 50)
            df_pca_role$Macro_State <- factor(paste("State", km_final$cluster, target_role, target_cond, sep="_"))
          } else {
            cluster_data <- pca_res$x[, 1:min(4, ncol(pca_res$x))]
            set.seed(42)
            gmm_res <- Mclust(cluster_data)
            optimal_k <- gmm_res$G
            df_pca_role$Macro_State <- factor(paste("State", gmm_res$classification, target_role, target_cond, sep="_"))
          }
          
          state_counts <- df_pca_role %>% group_by(Macro_State) %>% summarise(n=n(), .groups="drop")
          threshold <- max(10, nrow(df_pca_role) * 0.05)
          small_states <- as.character(state_counts %>% filter(n < threshold) %>% pull(Macro_State))
          
          if(length(small_states) > 0 && length(unique(df_pca_role$Macro_State)) > 1) {
            centroids <- df_pca_role %>% group_by(Macro_State) %>% summarise(across(all_of(valid_dist_cols), mean)) %>% column_to_rownames("Macro_State")
            for (ss in small_states) {
              valid_targets <- setdiff(rownames(centroids), small_states)
              if(length(valid_targets) == 0) break
              if(length(valid_targets) == 1) nearest <- valid_targets else nearest <- names(which.min(as.matrix(dist(centroids))[ss, valid_targets]))
              df_pca_role <- df_pca_role %>% mutate(Macro_State = as.character(Macro_State), Macro_State = ifelse(Macro_State == ss, nearest, Macro_State))
            }
            df_pca_role$Macro_State <- factor(df_pca_role$Macro_State)
            optimal_k <- length(unique(df_pca_role$Macro_State))
          }

          df_pca_role$Display_State <- factor(str_remove(as.character(df_pca_role$Macro_State), paste0("_", target_role, "_", target_cond)))
          
          df_pca_role <- df_pca_role %>% mutate(Test_Group = paste(grp, Complex_Type, sep="@"))
          tbl_var <- table(df_pca_role$Test_Group, df_pca_role$Macro_State)
          if(nrow(tbl_var) > 1 && ncol(tbl_var) > 1) {
            pval_var <- tryCatch(fisher.test(tbl_var, simulate.p.value=TRUE, B=5000)$p.value, error=function(e) NA)
            write_csv(tibble(Target_Role=target_role, Cond_State=target_cond, Global_p_value=pval_var) %>% mutate(Significance = case_when(Global_p_value < 0.001 ~ "***", Global_p_value < 0.01 ~ "**", Global_p_value < 0.05 ~ "*", TRUE ~ "ns")), file.path(family_out_dir, sprintf("Stats_P6_State_Variant_GLOBAL_%s_%s.csv", target_role, target_cond)))
            
            pw_var <- run_safe_pairwise_fisher(as.matrix(tbl_var))
            if(!is.null(pw_var)) {
               write_stats_csvs(format_pairwise_stats(pw_var, is_directional=FALSE), family_out_dir, sprintf("Stats_P6_State_Variant_%s_%s", target_role, target_cond), is_directional=FALSE)
            }
          }
          
          state_colors <- get_distinct_colors(optimal_k)
          names(state_colors) <- levels(df_pca_role$Display_State)
          
          plot_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], State = df_pca_role$Display_State)
          hulls <- plot_df %>% group_by(State) %>% slice(chull(PC1, PC2)) %>% ungroup()
          
          p_pca <- ggplot(plot_df, aes(x = PC1, y = PC2, color = State, fill = State)) + geom_point(alpha = 0.7, size = 2) + geom_polygon(data = hulls, alpha = 0.2, linewidth = 0) + scale_color_manual(values = state_colors) + scale_fill_manual(values = state_colors) + labs(title = sprintf("%s Meta-States (%s | %s)", grp, target_role, target_cond)) + PLOT_THEME
          suppressWarnings(ggsave(file.path(family_out_dir, sprintf("Plot_P6_State_Clusters_PCA_%s_%s.pdf", target_role, target_cond)), p_pca, width = 10, height = 7))
          
          AXIS_STYLE <- list(titlefont = list(size = 18, color = "black"), tickfont = list(size = 14, color = "black"))
          p_3d <- plot_ly(data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2], PC3=pca_res$x[,3], State=df_pca_role$Display_State), x=~PC1, y=~PC2, z=~PC3, color=~State, colors=state_colors, type='scatter3d', mode='markers', marker=list(size=4, opacity=0.8, line=list(width=0))) %>%
            plotly::layout(title = list(text = sprintf("<b>%s Phase Space (%s | %s)</b>", grp, target_role, target_cond), font = list(size = 24, color = "black")), scene = list(xaxis = c(list(title = 'PC1'), AXIS_STYLE), yaxis = c(list(title = 'PC2'), AXIS_STYLE), zaxis = c(list(title = 'PC3'), AXIS_STYLE)))
          suppressWarnings(htmlwidgets::saveWidget(p_3d, file.path(family_out_dir, sprintf("Plot_P6_Interactive_3D_Space_%s_%s.html", target_role, target_cond))))
          
          p_variant_state <- df_pca_role %>% group_by(Display_State, Complex_Type) %>% summarise(n = n(), .groups="drop") %>% ggplot(aes(x = Display_State, y = n, fill = Complex_Type)) + geom_bar(stat="identity", position="stack", color="black") + PLOT_THEME + labs(title = sprintf("Variant Composition of States (%s | %s)", target_role, target_cond), x="Meta-State", y = "Count", fill = "Complex Type") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
          suppressWarnings(ggsave(file.path(family_out_dir, sprintf("Plot_P6_State_Variant_Composition_%s_%s.pdf", target_role, target_cond)), p_variant_state, width=10, height=6))
          
          for(cv in names(TARGET_CATEGORIES)) df_pca_role[[paste0("Binary_", cv)]] <- factor(if_else(df_pca_role[[cv]]==TARGET_CATEGORIES[[cv]], TARGET_CATEGORIES[[cv]], "Other"), levels=c(TARGET_CATEGORIES[[cv]], "Other"))
          
          plot_feat <- function(v_name, o_var) {
            df_pca_role %>% group_by(Display_State, !!sym(v_name)) %>% summarise(n = n(), .groups="drop") %>% group_by(Display_State) %>% mutate(Percent=n/sum(n)*100) %>%
              ggplot(aes(x=Display_State, y=Percent, fill=!!sym(v_name))) + geom_bar(stat="identity", color="black", width=0.6) + scale_fill_manual(values=setNames(c("#2171b5", "#cccccc"), c(TARGET_CATEGORIES[[o_var]], "Other"))) + labs(title=sprintf("%s Identity", o_var), x="Meta-State", y="Proportion (%)", fill="Feature") + PLOT_THEME
          }
          
          p_kincore <- (plot_feat("Binary_State","State") | plot_feat("Binary_C_Helix","C_Helix")) / (plot_feat("Binary_R_Spine","R_Spine") | plot_feat("Binary_Spatial","Spatial")) + plot_annotation(title = sprintf("Biological Signatures (%s | %s)", target_role, target_cond), theme = theme(plot.title = element_text(size = 18, face = "bold")))
          suppressWarnings(ggsave(file.path(family_out_dir, sprintf("Plot_P7_MacroState_Signatures_%s_%s.pdf", target_role, target_cond)), p_kincore, width=16, height=12))
          
          df_pca_role <- df_pca_role %>% select(-Plot_Network, -Display_State, -Test_Group)
          df_pca_master_out <- bind_rows(df_pca_master_out, df_pca_role)
        }
      }
      
      if(nrow(df_pca_master_out) > 0) {
          write_csv(df_pca_master_out, file.path(family_out_dir, "Phase7_Complete_Structural_Metadata.csv"))
          ultimate_pca_master_out <- bind_rows(ultimate_pca_master_out, df_pca_master_out)
      }
    }
  }
}

if(nrow(ultimate_pca_master_out) > 0) {
  ultimate_pca_master_out <- ultimate_pca_master_out %>% distinct(Simulation_ID, Chain, .keep_all = TRUE)
  write_csv(ultimate_pca_master_out, file.path(OUT_DIR, "Phase7_Complete_Structural_Metadata.csv"))
  cat(sprintf("\n[✓] Master Phase 7 structural metadata generated. Ready for Phase 8.\n"))
}
