suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(scales)
})

options(bitmapType = "cairo")

# =========================================================
# paths
# =========================================================
BASE_DIR <- "/rds/general/project/hda_25-26/live/TDS/TDS_Group1"

REFIT_SPLIT_DIR <- file.path(
  BASE_DIR, "analysis", "mediation", "outputs", "formal_mediation_refit_split"
)

SHORTLIST_DIR <- file.path(
  BASE_DIR, "analysis", "mediation", "outputs", "final_shortlists"
)

OUT_DIR <- file.path(
  BASE_DIR, "analysis", "mediation", "outputs", "final_mediation_outputs"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

ANALYSES_TO_KEEP <- c("main", "female")

# =========================================================
# helper: safely read csv
# =========================================================
safe_read_csv <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  read_csv(path, show_col_types = FALSE)
}

# =========================================================
# helper: find chunk result files for one analysis
# =========================================================
find_result_files <- function(refit_split_dir, analysis_name) {
  pattern <- paste0("^formal_mediation_results_", analysis_name, "(_chunk[0-9]+)?\\.csv$")
  files <- list.files(
    refit_split_dir,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  sort(files)
}

# =========================================================
# helper: combine chunk result files
# =========================================================
combine_analysis_results <- function(refit_split_dir, analysis_name) {
  files <- find_result_files(refit_split_dir, analysis_name)
  
  if (length(files) == 0) {
    warning("No result files found for analysis: ", analysis_name)
    return(NULL)
  }
  
  df <- bind_rows(lapply(files, safe_read_csv)) %>%
    mutate(
      analysis = analysis_name,
      pathway = paste(exposure_term, mediator, sep = " -> ")
    ) %>%
    distinct()
  
  req_cols <- c(
    "analysis", "subgroup_value", "exposure_var", "exposure_term", "mediator",
    "n_complete",
    "a_estimate", "a_se", "a_p",
    "b_estimate", "b_se", "b_p",
    "total_effect_c", "total_se", "total_p",
    "direct_effect_cprime", "direct_se", "direct_p",
    "indirect_effect_ab", "indirect_direction", "note",
    "indirect_boot_mean", "indirect_ci_low", "indirect_ci_high",
    "direct_boot_mean", "direct_ci_low", "direct_ci_high",
    "total_boot_mean", "total_ci_low", "total_ci_high",
    "pathway"
  )
  
  missing_cols <- setdiff(req_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in ", analysis_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df %>%
    arrange(exposure_term, mediator)
}

# =========================================================
# formatting helpers
# =========================================================
fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), NA_character_, formatC(x, digits = digits, format = "fg", flag = "#"))
}

fmt_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

fmt_ci <- function(est, low, high, digits = 4) {
  ifelse(
    is.na(est) | is.na(low) | is.na(high),
    NA_character_,
    paste0(
      fmt_num(est, digits), " [",
      fmt_num(low, digits), ", ",
      fmt_num(high, digits), "]"
    )
  )
}

# =========================================================
# make publication-style table dataframe
# =========================================================
make_display_table <- function(df) {
  df %>%
    mutate(
      Direction = case_when(
        indirect_direction %in% c("positive", "negative") ~ indirect_direction,
        !is.na(indirect_boot_mean) & indirect_boot_mean > 0 ~ "positive",
        !is.na(indirect_boot_mean) & indirect_boot_mean < 0 ~ "negative",
        TRUE ~ NA_character_
      )
    ) %>%
    transmute(
      Exposure = exposure_term,
      Mediator = mediator,
      `N` = n_complete,
      `a path\n(est, p)` = paste0(fmt_num(a_estimate), ", ", fmt_p(a_p)),
      `b path\n(est, p)` = paste0(fmt_num(b_estimate), ", ", fmt_p(b_p)),
      `Indirect effect\nboot mean [95% CI]` = fmt_ci(indirect_boot_mean, indirect_ci_low, indirect_ci_high),
      `Direct effect\nboot mean [95% CI]` = fmt_ci(direct_boot_mean, direct_ci_low, direct_ci_high),
      `Total effect\nboot mean [95% CI]` = fmt_ci(total_boot_mean, total_ci_low, total_ci_high),
      `Direction` = Direction
    )
}

# =========================================================
# save PNG table
# =========================================================
save_table_png <- function(table_df, title_text, out_path, width = 22, height = 4.2, base_size = 11) {
  table_theme <- ttheme_minimal(
    base_size = base_size,
    core = list(
      fg_params = list(hjust = 0.5, x = 0.5, fontsize = base_size),
      padding = unit(c(6, 6), "mm")
    ),
    colhead = list(
      fg_params = list(fontface = 2, hjust = 0.5, x = 0.5, fontsize = base_size),
      padding = unit(c(7, 7), "mm")
    )
  )
  
  tg <- tableGrob(table_df, rows = NULL, theme = table_theme)
  
  tg$widths <- unit(
    c(2.3, 2.5, 0.8, 1.8, 1.8, 3.4, 3.4, 3.4, 1.2),
    "in"
  )
  
  png(out_path, width = width, height = height, units = "in", res = 300, type = "cairo")
  grid.newpage()
  grid.draw(
    arrangeGrob(
      textGrob(
        title_text,
        gp = gpar(fontsize = 18, fontface = "bold"),
        x = 0.5,
        hjust = 0.5
      ),
      tg,
      ncol = 1,
      heights = unit.c(unit(0.7, "in"), unit(1, "npc") - unit(0.7, "in"))
    )
  )
  dev.off()
}

# =========================================================
# forest plot for indirect effects
# =========================================================
make_forest_plot <- function(df, title_text, out_path) {
  plot_df <- df %>%
    mutate(
      pathway = paste(exposure_term, mediator, sep = " -> "),
      pathway_label = stringr::str_wrap(pathway, width = 28),
      pathway_label = factor(pathway_label, levels = rev(unique(pathway_label)))
    )
  
  p <- ggplot(
    plot_df,
    aes(x = indirect_boot_mean, y = pathway_label)
  ) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.7) +
    geom_errorbarh(
      aes(xmin = indirect_ci_low, xmax = indirect_ci_high),
      height = 0.18,
      linewidth = 0.9
    ) +
    geom_point(size = 3) +
    labs(
      title = title_text,
      x = "Indirect effect (bootstrap mean, 95% CI)",
      y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.text.y = element_text(size = 11),
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 20, 15, 30)
    )
  
  ggsave(out_path, p, width = 11, height = 6, dpi = 300, bg = "white")
}

# =========================================================
# pathway counts summary
# =========================================================
make_counts_summary <- function(main_df, female_df, shortlist_dir) {
  male_n <- 0L
  
  male_summary_path <- file.path(shortlist_dir, "summary_shortlist_male.csv")
  if (file.exists(male_summary_path)) {
    male_summary <- safe_read_csv(male_summary_path)
    if ("n_strict_shortlist" %in% names(male_summary)) {
      male_n <- as.integer(male_summary$n_strict_shortlist[[1]])
    }
  }
  
  tibble(
    analysis = c("main", "female", "male"),
    n_pathways = c(
      if (is.null(main_df)) 0L else nrow(main_df),
      if (is.null(female_df)) 0L else nrow(female_df),
      male_n
    )
  )
}

make_counts_plot <- function(count_df, out_path) {
  p <- ggplot(count_df, aes(x = analysis, y = n_pathways)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n_pathways), vjust = -0.35, size = 5) +
    labs(
      title = "Shortlisted mediation pathways by analysis",
      x = NULL,
      y = "Number of pathways"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      panel.grid.minor = element_blank(),
      plot.margin = margin(15, 20, 15, 20)
    ) +
    expand_limits(y = max(count_df$n_pathways) + 1)
  
  ggsave(out_path, p, width = 8, height = 5.5, dpi = 300, bg = "white")
}

# =========================================================
# layered mediation pathway figure
# =========================================================
draw_layered_panel <- function(df, panel_title, show_effect = FALSE) {
  if (is.null(df) || nrow(df) == 0) {
    plot.new()
    title(main = panel_title, cex.main = 1.4, font.main = 2)
    text(0.5, 0.5, "No pathways available", cex = 1.1)
    return(invisible(NULL))
  }
  
  exposures <- unique(df$exposure_term)
  mediators <- unique(df$mediator)
  outcome <- "cvd_incident"
  
  exp_y <- seq(0.68, 0.22, length.out = length(exposures))
  med_y <- seq(0.72, 0.18, length.out = length(mediators))
  out_y <- 0.5
  
  exp_pos <- data.frame(name = exposures, x = 0.18, y = exp_y, stringsAsFactors = FALSE)
  med_pos <- data.frame(name = mediators, x = 0.52, y = med_y, stringsAsFactors = FALSE)
  out_pos <- data.frame(name = outcome, x = 0.84, y = out_y, stringsAsFactors = FALSE)
  
  clean_label <- function(x) {
    stringr::str_wrap(gsub("_", " ", x), width = 14)
  }
  
  get_col <- function(direction) {
    ifelse(direction %in% "negative", "firebrick3", "darkgreen")
  }
  
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title(main = panel_title, cex.main = 1.5, font.main = 2)
  
  text(0.18, 0.985, "Exposure", font = 2, cex = 1.1)
  text(0.52, 0.985, "Mediator", font = 2, cex = 1.1)
  text(0.84, 0.985, "Outcome", font = 2, cex = 1.1)
  
  for (i in seq_len(nrow(df))) {
    ex <- df$exposure_term[i]
    md <- df$mediator[i]
    dir <- ifelse(
      !is.na(df$indirect_direction[i]),
      df$indirect_direction[i],
      ifelse(df$indirect_boot_mean[i] < 0, "negative", "positive")
    )
    
    ex_row <- exp_pos[exp_pos$name == ex, , drop = FALSE]
    md_row <- med_pos[med_pos$name == md, , drop = FALSE]
    
    arrows(
      x0 = ex_row$x + 0.07, y0 = ex_row$y,
      x1 = md_row$x - 0.07, y1 = md_row$y,
      length = 0.08, lwd = 1.6, col = get_col(dir)
    )
    
    arrows(
      x0 = md_row$x + 0.08, y0 = md_row$y,
      x1 = out_pos$x - 0.08, y1 = out_pos$y,
      length = 0.08, lwd = 1.6, col = get_col(dir)
    )
    
    if (show_effect) {
      label_x <- (md_row$x + out_pos$x) / 2
      label_y <- (md_row$y + out_pos$y) / 2 + 0.03
      text(
        label_x, label_y,
        labels = sprintf("%.2e", df$indirect_boot_mean[i]),
        cex = 0.75, col = get_col(dir)
      )
    }
  }
  
  draw_nodes <- function(df_nodes, bg = "lightblue", radius = 0.058) {
    for (i in seq_len(nrow(df_nodes))) {
      symbols(
        x = df_nodes$x[i], y = df_nodes$y[i],
        circles = radius, inches = FALSE, add = TRUE,
        bg = bg, fg = "grey35"
      )
      text(
        df_nodes$x[i], df_nodes$y[i],
        labels = clean_label(df_nodes$name[i]),
        cex = 0.64
      )
    }
  }
  
  draw_nodes(exp_pos, bg = "#a6cee3", radius = 0.055)
  draw_nodes(med_pos, bg = "#a6cee3", radius = 0.050)
  draw_nodes(out_pos, bg = "#f1d9a7", radius = 0.055)
  
  legend(
    "bottomleft",
    legend = c("positive indirect direction", "negative indirect direction"),
    col = c("darkgreen", "firebrick3"),
    lwd = 2,
    bty = "n",
    cex = 0.95
  )
}

save_layered_single <- function(df, panel_title, out_path, show_effect = FALSE) {
  png(out_path, width = 10, height = 7.5, units = "in", res = 300, type = "cairo")
  par(mfrow = c(1, 1), mar = c(2.5, 2.5, 4, 2.5))
  draw_layered_panel(df, panel_title, show_effect = show_effect)
  dev.off()
}

# =========================================================
# main
# =========================================================
main_df <- combine_analysis_results(REFIT_SPLIT_DIR, "main")
female_df <- combine_analysis_results(REFIT_SPLIT_DIR, "female")

# save merged raw csv
write_csv(main_df, file.path(OUT_DIR, "final_mediation_results_main.csv"))
write_csv(female_df, file.path(OUT_DIR, "final_mediation_results_female.csv"))

# display tables
main_table <- make_display_table(main_df)
female_table <- make_display_table(female_df)

write_csv(main_table, file.path(OUT_DIR, "display_table_main.csv"))
write_csv(female_table, file.path(OUT_DIR, "display_table_female.csv"))

save_table_png(
  main_table,
  "Main analysis: shortlisted exposure–mediator–outcome pathways",
  file.path(OUT_DIR, "mediation_main_table.png"),
  width = 22,
  height = 4.8 + 0.45 * nrow(main_table),
  base_size = 11
)

save_table_png(
  female_table,
  "Female-stratified analysis: shortlisted exposure–mediator–outcome pathways",
  file.path(OUT_DIR, "mediation_female_table.png"),
  width = 22,
  height = 4.8 + 0.45 * nrow(female_table),
  base_size = 11
)


# forest plots
make_forest_plot(
  main_df,
  "Main analysis: indirect effect estimates",
  file.path(OUT_DIR, "forest_indirect_main.png")
)

make_forest_plot(
  female_df,
  "Female-stratified analysis: indirect effect estimates",
  file.path(OUT_DIR, "forest_indirect_female.png")
)

# counts summary
count_df <- make_counts_summary(main_df, female_df, SHORTLIST_DIR)
write_csv(count_df, file.path(OUT_DIR, "pathway_counts_summary.csv"))

make_counts_plot(
  count_df,
  file.path(OUT_DIR, "pathway_counts_summary.png")
)

# main only
save_layered_single(
  main_df,
  "Main analysis: layered mediation pathway figure",
  file.path(OUT_DIR, "layered_main.png"),
  show_effect = FALSE
)

save_layered_single(
  main_df,
  "Main analysis: layered mediation pathway figure",
  file.path(OUT_DIR, "layered_main_with_effects.png"),
  show_effect = TRUE
)

# female only
save_layered_single(
  female_df,
  "Female-stratified analysis: layered mediation pathway figure",
  file.path(OUT_DIR, "layered_female.png"),
  show_effect = FALSE
)

save_layered_single(
  female_df,
  "Female-stratified analysis: layered mediation pathway figure",
  file.path(OUT_DIR, "layered_female_with_effects.png"),
  show_effect = TRUE
)

cat("All mediation outputs saved to:\n", OUT_DIR, "\n")

