###############################################################################
# Load packages
load_packages <- function() {
  cat("Loading packages", sep = "\n")
  # Package dependencies
  pkgs <- c(
    "data.table",
    "tidyverse",
    "here",
    "ggplot2",
    "scales",
    "RColorBrewer",
    "ggpubr",
    "glmmTMB",
    "ggrastr",
    "patchwork",
    "readxl"
  )
  invisible(lapply(pkgs, library, character.only = TRUE))
}
suppressMessages(load_packages())


###############################################################################
outdir <- here("outputs")
cat("Saving output data and figures to: ", outdir, sep = "\n")
dir.create(outdir, showWarnings = TRUE)
###############################################################################


#### Factor order for variables ====
# Top level cell type ordering
knn_top_ord <- c(
  "TNKILC", "B", "Plasma", "Mast", "Myeloid", "Strom", "Epi"
)

# Mid level cell type ordering
knn_mid_ord <- c(
  "Tcd4", "Tcd8", "Tplzf", "NK", "ILC", "B", "Plasma", "Mast",
  "Myeloid", "Endo", "Peri", "Fibro",
  "SmoothMuscle", "Schwann", "Epi"
)

knn_group_ord <- c(
  "Immune", "Endothelial", "Pericyte", "Fibroblast", "SmoothMuscle", "Schwann", "Epithelial"
)

# Pathology annotations
path_ord <- c(
  "tumor", "tumor_invasive_margin", "tumor_luminal_margin",
  "non_neoplastic_mucosa", "non_neoplastic_other",
  "lymphoid_structure"
)

#### Data files ====
cat("Loading data files", sep = "\n")
# Single cell annotations
annots <- fread(here("data/GSE_files/single_cell_annots_obs.csv"))

# Patient metadata
patient_metadata <- fread(here("data/GSE_files/metatable.csv"))

# Celltype map
knn_top_map <- annots %>%
  select(KNN_Top, KNN_Mid, KNN_celltype_v2, KNN_Grouped) %>%
  distinct() %>%
  arrange(KNN_Top, KNN_Mid, KNN_celltype_v2) %>%
  mutate(
    KNN_Top = factor(KNN_Top, levels = knn_top_ord),
    KNN_Mid = factor(KNN_Mid, levels = knn_mid_ord)
  )

#### Patient and MMR status sample order ====
# MMR status sample order
sample_ord <- patient_metadata %>%
  select(File, MMRstatus) %>%
  arrange(MMRstatus, File) %>%
  pull(File)

patient_ord <- patient_metadata %>%
  select(PatientID, MMRstatus) %>%
  distinct() %>%
  arrange(MMRstatus, PatientID) %>%
  pull(PatientID)

#### Color maps ====
# Top level celltype annotations
cmap_top <- c(
  "TNKILC" = "#ff0000",
  "B" = "#0011ff",
  "Plasma" = "#66b5ff",
  "Mast" = "#f0e928",
  "Myeloid" = "#ffaf24",
  "Strom" = "#09ceca",
  "Epi" = "#ae00ff"
)

# Pathology annotations
cmap_path_annot <- c(
  "tumor" = "#8A9FD1",
  "tumor_invasive_margin" = "#1C75BC",
  "tumor_luminal_margin" = "#372461",
  "non_neoplastic_mucosa" = "#FAD14A",
  "non_neoplastic_other" = "#D3D3D3",
  "lymphoid_structure" = "#E68428"
)

# MMR status
cmap_MMR <- c(
  "MMRd" = "#BE1E2D",
  "MMRp" = "#1C75BC"
)


cmap_tessera_annot <- c(
  "stromal_tile" = "#18527e",
  "epithelial_tile" = "#a97f2f",
  "non_tumor" = "#961f1f"
)

cmap_band <- c(
  "thin" = "#427aa1",
  "thick" = "#064789"
)

cmap_CXCL <- c(
  "CXCL_pos" = "#660d20",
  "CXCL_neg" = "#898e9f"
)


###############################################################################
# Plotting aesthetics
set_figure_size <- function(
  width,
  height,
  res = 300
) {
  options(
    repr.plot.width = width,
    repr.plot.height = height,
    repr.plot.res = res
  )
}

# Default ggplot2 theme
theme_set(
  theme_bw() +
    theme(
      text = element_text(
        family = "Helvetica",
        color = "black",
        size = 8
      ),
      title = element_text(
        size = 8
      ),
      # Axis:
      axis.title = element_text(
        size = 8
      ),
      axis.text = element_text(
        size = 8,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.25
      ),
      # Plot:
      plot.background = element_blank(),
      # Panel:
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.25
      ),
      # Legend:
      legend.text = element_text(
        size = 8
      ),
      legend.title = element_text(
        size = 8
      ),
      legend.key.size = unit(
        0.4,
        "cm"
      ),
      legend.background = element_blank(),
      legend.key = element_blank(),
      # Facets:
      strip.background = element_blank(),
      strip.text = element_text(
        size = 8
      )
    )
)

# Vertically orient X axis
vertical_x_text <- function(...) {
  vertical_theme <- theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    validate = TRUE,
    ...
  )
  return(vertical_theme)
}

# Remove everything from x axis
strip_x_axis <- function(...) {
  strip_x_theme <- theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    validate = TRUE,
    ...
  )
  return(strip_x_theme)
}

# Remove everything from y axis
strip_y_axis <- function(...) {
  strip_y_theme <- theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    validate = TRUE,
    ...
  )
  return(strip_y_theme)
}


# Fill RdBu color map for heatmaps
scale_rdbu_fill <- function(limits,
                            ...) {
  # Check limits input to ensure there is a minimum and a maximum value
  stopifnot(length(limits) == 3)
  hmap.lim.min <- limits[1]
  hmap.lim.mid <- limits[2]
  hmap.lim.max <- limits[3]

  scale_fill_gradientn(
    colors = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    values = scales::rescale(c(
      hmap.lim.min,
      hmap.lim.mid,
      hmap.lim.max
    )),
    limits = c(hmap.lim.min, hmap.lim.max),
    oob = scales::squish,
    na.value = "lightgrey",
    ...
  )
}

###############################################################################
`%out%` <- function(a, b) {
  !a %in% b
}


get_heatmap_hclust <- function(
  df_input,
  x_col,
  y_col,
  fill_col,
  na_fill = NA,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  dist_method = "euclidean",
  hclust_method = "complete"
) {
  # Symbols
  x_sym <- rlang::sym(x_col)
  y_sym <- rlang::sym(y_col)
  fill_sym <- rlang::sym(fill_col)

  # Subset relevant columns
  df_input <- df_input %>%
    select(!!x_sym, !!y_sym, !!fill_sym) %>%
    ungroup()

  get_order <- function(mat, axis = "row") {
    dist_mat <- stats::dist(mat, method = dist_method)
    stats::hclust(dist_mat, method = hclust_method)$order
  }

  result <- list()

  if (isTRUE(cluster_rows)) {
    row_mat <- df_input %>%
      pivot_wider(
        names_from = !!x_sym,
        values_from = !!fill_sym,
        values_fill = na_fill
      )

    row_names <- row_mat[[1]]
    mat <- as.matrix(row_mat[, -1])
    result$y_order <- row_names[get_order(mat, "row")]
  }

  if (isTRUE(cluster_columns)) {
    col_mat <- df_input %>%
      pivot_wider(
        names_from = !!y_sym,
        values_from = !!fill_sym,
        values_fill = na_fill
      )

    col_names <- col_mat[[1]]
    mat <- as.matrix(col_mat[, -1])
    result$x_order <- col_names[get_order(mat, "col")]
  }

  # If only one clustering is requested, return vector instead of list
  if (cluster_rows && !cluster_columns) {
    return(result$y_order)
  }
  if (cluster_columns && !cluster_rows) {
    return(result$x_order)
  }
  return(result)
}


###############################################################################
rotate_xy <- function(
  df,
  x_col = "x", y_col = "y",
  angle_deg
) {
  # Convert degrees to radians
  theta <- angle_deg * pi / 180

  # Rotation matrix
  R <- matrix(
    c(
      cos(theta), -sin(theta),
      sin(theta), cos(theta)
    ),
    nrow = 2, byrow = TRUE
  )

  # Extract coordinates
  coords <- df %>%
    select(!!as.symbol(x_col), !!as.symbol(y_col)) %>%
    as.matrix()

  # Rotate
  rotated <- coords %*% t(R) %>%
    as.data.frame() %>%
    rename(x_rotated = V1, y_rotated = V2)

  # Replace in dataframe
  df <- bind_cols(df, rotated)

  return(df)
}


add_scale_bar <- function(
  p,
  scale_len = 1000,
  scale_width = 0.5,
  text_size = 1.5
) {
  # Get the original data
  df_original <- p$data
  aes_mapping <- p$mapping

  # Get the x and y variable names as strings
  x_var <- rlang::as_label(aes_mapping$x)
  y_var <- rlang::as_label(aes_mapping$y)

  # Now you can extract the actual x and y values from the data
  df_original <- p$data
  x_vals <- df_original[[x_var]]
  y_vals <- df_original[[y_var]]

  # Get min/max
  x_min <- min(x_vals, na.rm = TRUE)
  x_max <- max(x_vals, na.rm = TRUE)
  y_min <- min(y_vals, na.rm = TRUE)
  y_max <- max(y_vals, na.rm = TRUE)
  y_expand <- y_max - y_min * 0.001

  # Add scale bar and text
  p <- p +
    annotate(
      "segment",
      x = x_min, xend = x_min + scale_len,
      y = y_min, yend = y_min,
      linewidth = scale_width
    ) +
    # Dynamically expand limits to ensure the scale bar fits and make equal
    coord_equal(
      xlim = c(x_min, x_max),
      ylim = c(y_min - y_expand, y_max)
    )

  if (text_size > 0) {
    p <- p +
      annotate(
        "text",
        x = x_min + scale_len / 2, y = y_min,
        label = paste0(scale_len, " \u03BCm"),
        vjust = 1.5,
        size = 1.5
      )
  }

  return(p)
}


plot_glmm_barplot <- function(
  glmm_res
) {
  plt <- glmm_res %>%
    filter(!grepl("intercept", rowname)) %>%
    left_join(
      knn_top_map,
      by = c("rowname" = "KNN_celltype_v2")
    ) %>%
    mutate(
      rowname = fct_reorder(rowname, log2fc),
      padj = p.adjust(Pr.z, method = "fdr")
    ) %>%
    ggplot(mapping = aes(x = rowname, y = log2fc, fill = KNN_Top)) +
    geom_bar(stat = "identity") +
    geom_text(
      data = . %>% filter(padj < 0.05),
      mapping = aes(y = log2fc + ifelse(log2fc > 0, 0.15, -0.25)),
      label = "*", size = 3
    ) +
    scale_fill_manual(values = cmap_top) +
    labs(
      x = NULL,
      y = "LogFC of GLMM effect size",
      fill = NULL
    ) +
    vertical_x_text()

  return(plt)
}
