# analysis/table1/render_table1_png_base.R
# Render Table1 CSV outputs to publication-style PNGs using base R only
# No gt/htmlwidgets required.

suppressWarnings({
  if (!requireNamespace("grid", quietly = TRUE)) stop("Base package 'grid' missing (unexpected).")
})

library(grid)

IN_DIR  <- file.path("analysis", "table1", "output")
OUT_DIR <- file.path("analysis", "table1", "output", "output_rendered_png")  # same as your old folder
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 1) DELETE old PNGs so only updated outputs remain
old_pngs <- list.files(OUT_DIR, pattern = "\\.png$", full.names = TRUE)
if (length(old_pngs) > 0) file.remove(old_pngs)

# Utility: safe numeric for sorting page names
pad_int <- function(x, w = 2) sprintf(paste0("%0", w, "d"), x)

# ---- Table drawing parameters (tweak if needed) ----
PAGE_W_IN  <- 11     # width in inches
PAGE_H_IN  <- 8.5    # height in inches
DPI        <- 220

# layout
LEFT_MARGIN  <- 0.02
RIGHT_MARGIN <- 0.02
TOP_MARGIN   <- 0.04
BOTTOM_MARGIN<- 0.04

# row/col sizing
BASE_FONTSIZE <- 9
HEADER_FONTSIZE <- 9.5

# Try Times; if not available, grid will fallback automatically
BASE_FONTFAMILY <- "Times"

# Decide how many rows per page (dynamic based on height/font)
rows_per_page <- function(n_rows) {
  # conservative default; you can increase if you want denser tables
  32
}

# Column width weights: make "Variable"/"level" wider
# Works with your CSV schema: variable, level, Overall, group1, group2, p_value, (optional missing col)
width_weights_for <- function(cols) {
  w <- rep(1, length(cols))
  names(w) <- cols
  if ("variable" %in% cols) w[cols == "variable"] <- 2.2
  if ("level" %in% cols)    w[cols == "level"]    <- 2.4
  # common value columns
  for (nm in c("Overall","No incident CVD","Incident CVD","Female","Male","p_value","Missing n (%) (Overall)")) {
    if (nm %in% cols) w[cols == nm] <- 1.2
  }
  # allow long label columns to be slightly wider
  w
}

# Word wrap helper (simple): inserts \n if too long
wrap_text <- function(x, width = 28) {
  x <- ifelse(is.na(x), "", as.character(x))
  vapply(x, function(s) {
    if (nchar(s) <= width) return(s)
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1))
}

draw_table_page <- function(df_page, title, out_png, page_i, page_n) {
  # Prepare display df
  df_disp <- df_page
  
  # Clean NAs
  df_disp[is.na(df_disp)] <- ""
  
  # Wrap long text in variable/level
  if ("variable" %in% names(df_disp)) df_disp$variable <- wrap_text(df_disp$variable, 22)
  if ("level" %in% names(df_disp))    df_disp$level    <- wrap_text(df_disp$level, 28)
  
  cols <- names(df_disp)
  nrowp <- nrow(df_disp)
  ncolp <- ncol(df_disp)
  
  # Column widths (normalized)
  ww <- width_weights_for(cols)
  col_w <- ww / sum(ww)
  
  # Row heights (header + body)
  header_h <- 0.06
  body_h   <- (1 - TOP_MARGIN - BOTTOM_MARGIN - header_h) / max(1, nrowp)
  
  # Start device
  png(out_png, width = PAGE_W_IN, height = PAGE_H_IN, units = "in", res = DPI)
  on.exit(dev.off(), add = TRUE)
  
  grid.newpage()
  
  # Title
  grid.text(
    label = paste0(title, if (page_n > 1) paste0(" (page ", page_i, "/", page_n, ")") else ""),
    x = LEFT_MARGIN, y = 1 - TOP_MARGIN/2,
    just = c("left", "top"),
    gp = gpar(fontsize = 14, fontface = "bold", fontfamily = BASE_FONTFAMILY)
  )
  
  # Table viewport
  vp <- viewport(
    x = 0.5, y = 0.5,
    width = 1 - LEFT_MARGIN - RIGHT_MARGIN,
    height = 1 - TOP_MARGIN - BOTTOM_MARGIN,
    just = c("center", "center")
  )
  pushViewport(vp)
  
  # Compute x positions
  x_lefts <- cumsum(c(0, col_w[-length(col_w)]))
  x_cent  <- x_lefts + col_w/2
  
  # Header background
  grid.rect(
    x = 0.5, y = 1 - header_h/2,
    width = 1, height = header_h,
    gp = gpar(fill = "grey95", col = "grey80")
  )
  
  # Header text
  for (j in seq_len(ncolp)) {
    lab <- cols[j]
    # nicer names if your CSV uses internal names
    # your outputs already use "No incident CVD" etc, so keep as-is
    grid.text(
      label = lab,
      x = x_cent[j], y = 1 - header_h/2,
      gp = gpar(fontsize = HEADER_FONTSIZE, fontface = "bold", fontfamily = BASE_FONTFAMILY)
    )
  }
  
  # Body rows
  for (i in seq_len(nrowp)) {
    y_top <- 1 - header_h - (i-1)*body_h
    y_mid <- y_top - body_h/2
    
    # alternating row shading
    if (i %% 2 == 0) {
      grid.rect(
        x = 0.5, y = y_mid,
        width = 1, height = body_h,
        gp = gpar(fill = "white", col = NA)
      )
    }
    
    # horizontal line
    grid.lines(
      x = unit(c(0,1), "npc"), y = unit(rep(y_top,2), "npc"),
      gp = gpar(col = "grey85", lwd = 0.6)
    )
    
    for (j in seq_len(ncolp)) {
      val <- as.character(df_disp[i, j])
      
      # alignment: left for variable/level, center otherwise
      is_left <- cols[j] %in% c("variable","level")
      grid.text(
        label = val,
        x = if (is_left) x_lefts[j] + 0.01 else x_cent[j],
        y = y_mid,
        just = if (is_left) c("left","center") else c("center","center"),
        gp = gpar(fontsize = BASE_FONTSIZE, fontfamily = BASE_FONTFAMILY)
      )
    }
  }
  
  # bottom line
  grid.lines(
    x = unit(c(0,1), "npc"), y = unit(rep(0,2), "npc"),
    gp = gpar(col = "grey70", lwd = 0.8)
  )
  
  popViewport()
}

render_one_csv <- function(csv_path) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  
  # If user wants NA% column in BEFORE tables, it should already be in CSV after you update pipeline.
  # But if not present and file name indicates "before", we can compute missing column here as fallback.
  is_before <- grepl("_before", basename(csv_path))
  if (is_before && !("Missing n (%) (Overall)" %in% names(df))) {
    # try to compute missing overall from the raw 'variable' columns? not possible from summary csv
    # so we only warn; proper way is to add in Table1 generation step.
    message("WARN: ", basename(csv_path), " has no 'Missing n (%) (Overall)' column. ",
            "You asked to add NA% in BEFORE tables — this must be added in Table1 pipeline, not at render stage.")
  }
  
  title <- tools::file_path_sans_ext(basename(csv_path))
  out_prefix <- file.path(OUT_DIR, title)
  
  rpp <- rows_per_page(nrow(df))
  pages <- split(df, ceiling(seq_len(nrow(df)) / rpp))
  page_n <- length(pages)
  
  for (k in seq_along(pages)) {
    out_png <- paste0(out_prefix, "_page", pad_int(k, 2), ".png")
    draw_table_page(pages[[k]], title = title, out_png = out_png, page_i = k, page_n = page_n)
  }
}

# ---- Run ----
csv_files <- list.files(IN_DIR, pattern = "^table1_.*\\.csv$", full.names = TRUE)
if (length(csv_files) == 0) stop("No table1_*.csv found under: ", IN_DIR)

message("Found ", length(csv_files), " CSVs. Rendering to: ", OUT_DIR)
for (f in csv_files) {
  message("Rendering: ", basename(f))
  render_one_csv(f)
}
message("Done. PNGs in: ", OUT_DIR)




# ==============================
# EXTRA: export one LONG PNG per CSV (no pagination)
# ==============================
suppressPackageStartupMessages({
  library(grid)
  library(gridExtra)
})

CSV_DIR  <- file.path("analysis", "table1", "output")
PNG_DIR  <- file.path(CSV_DIR, "output_rendered_png")

# long image sizing (tune if needed)
DPI <- 220
W_IN <- 11.0         # width in inches
ROW_IN <- 0.22       # inches per row (0.18~0.28)
MIN_H <- 8.5
MAX_H <- 60

csv_files <- list.files(CSV_DIR, pattern = "\\.csv$", full.names = TRUE)
csv_files <- csv_files[grepl("^table1_", basename(csv_files))]

for (csv_path in csv_files) {
  df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) next
  
  # replace NA with blank so table looks clean
  df[is.na(df)] <- ""
  
  base <- sub("\\.csv$", "", basename(csv_path))
  out_long <- file.path(PNG_DIR, paste0(base, "_LONG.png"))
  
  # auto height based on number of rows
  h_in <- max(MIN_H, min(MAX_H, (nrow(df) + 6) * ROW_IN))
  
  png(out_long, width = W_IN, height = h_in, units = "in", res = DPI)
  grid.newpage()
  
  # a simple title
  grid.text(base, x = 0.02, y = 0.99, just = c("left", "top"),
            gp = gpar(fontsize = 14, fontface = "bold"))
  
  # table body
  tbl <- tableGrob(df, rows = NULL,
                   theme = ttheme_minimal(
                     base_size = 9,
                     core = list(padding = unit(c(2, 2), "mm")),
                     colhead = list(fg_params = list(fontface = "bold"))
                   ))
  
  # push table slightly down to leave space for title
  vp <- viewport(x = 0.5, y = 0.48, width = 0.98, height = 0.93)
  pushViewport(vp)
  grid.draw(tbl)
  upViewport()
  
  dev.off()
  message("Wrote LONG: ", out_long)
}