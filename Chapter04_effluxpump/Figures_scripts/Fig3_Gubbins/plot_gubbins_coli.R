
getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 04 Campy cmeABC efflux pump/WRITING/Figures/FIG_gubbins/")

###############################################################################
## Visualise Gubbins recombination in the cmeABCR region                     ##
## Aligned to phylogenetic tree with gene annotation, metadata & heatmap     ##
###############################################################################

# ── Libraries ────────────────────────────────────────────────────────────────
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(gggenes)
library(aplot)
library(patchwork)
library(cowplot)
library(RColorBrewer)


# ── INPUT FILES ──────────────────────────────────────────────────────────────
tree_file    <- "./coli_gubbins_tree.tre"
rec_gff_file <- "coli_gubbins_recombination.gff"
meta_file    <- "Metadata_coli.csv"
output_file  <- "c_coli_cmeABCR_recombination.pdf"


# ── GENE COORDINATES ────────────────────────────────────────────────────────
gene_df <- tibble::tribble(
  ~gene,   ~start,  ~end,    ~strand, ~direction,
  "cmeC",  373168,  374646,  -1,      -1,
  "cmeB",  374639,  377761,  -1,      -1,
  "cmeA",  377761,  378864,  -1,      -1,
  "cmeR",  378972,  379607,  -1,      -1
)

padding      <- 500
region_start <- min(gene_df$start) - padding
region_end   <- max(gene_df$end)   + padding


# ── OUTPUT SETTINGS ──────────────────────────────────────────────────────────
output_width  <- 16
output_height <- 10


# ── COLOUR PALETTES ─────────────────────────────────────────────────────────
rec_colours <- c("shared" = "#D73027", "unique" = "#4575B4")
gene_fills  <- c("cmeR" = "#FDB462", "cmeA" = "#80B1D3",
                 "cmeB" = "#FB8072", "cmeC" = "#B3DE69")

many_colour_palette <- colorRampPalette(
  c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
    "#ff7f00", "#a65628", "#f781bf", "#999999",
    "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
    "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
)


###############################################################################
## 1. READ AND PROCESS TREE                                                  ##
###############################################################################

tree   <- ape::read.tree(tree_file)
n_taxa <- length(tree$tip.label)

tree_linewidth <- case_when(
  n_taxa < 20  ~ 0.8,
  n_taxa < 50  ~ 0.5,
  n_taxa < 100 ~ 0.35,
  TRUE         ~ 0.25
)

tip_label_size <- case_when(
  n_taxa < 20  ~ 3.5,
  n_taxa < 50  ~ 2.5,
  n_taxa < 100 ~ 2.0,
  TRUE         ~ 1.5
)

p_tree <- ggtree(tree, size = tree_linewidth) +
  # geom_tiplab removed — no tip labels
  theme_tree2() +
  theme(
    legend.position = "none",
    axis.line.x     = element_line(linewidth = 0.25),
    plot.margin     = margin(5, 0, 5, 5)
  )

# Get taxa order from tree (important for aligning all panels)
taxa_order <- rev(get_taxa_name(p_tree))


###############################################################################
## 2. PARSE GUBBINS RECOMBINATION GFF                                        ##
###############################################################################

parse_gubbins_gff <- function(gff_file) {
  raw <- readLines(gff_file)
  raw <- raw[!grepl("^##", raw)]
  
  df <- read.table(text = raw, sep = "\t", header = FALSE,
                   stringsAsFactors = FALSE, fill = TRUE, quote = "")
  colnames(df) <- c("seq", "prog", "type", "start", "end",
                    "score", "strand", "frame", "attributes")
  
  df <- df %>%
    rowwise() %>%
    mutate(
      taxa_str = str_extract(attributes, "taxa=\\s*[^;]+"),
      taxa_str = str_remove(taxa_str, "taxa=\\s*"),
      taxa_str = str_trim(taxa_str),
      taxa     = list(str_split(taxa_str, "\\s+")[[1]])
    ) %>%
    ungroup() %>%
    mutate(
      n_taxa   = map_int(taxa, length),
      rec_type = if_else(n_taxa == 1, "unique", "shared")
    ) %>%
    select(start, end, taxa, n_taxa, rec_type)
  
  return(df)
}

rec_df <- parse_gubbins_gff(rec_gff_file) %>%
  mutate(
    start = as.numeric(start),
    end   = as.numeric(end)
  )

# Filter recombinations overlapping the cmeABCR region
rec_region <- rec_df %>%
  filter(end >= region_start & start <= region_end) %>%
  mutate(
    start = pmax(start, region_start),
    end   = pmin(end, region_end)
  )

# Unnest to one row per taxon
rec_long <- rec_region %>%
  unnest(taxa) %>%
  rename(label = taxa) %>%
  filter(!is.na(start), !is.na(end)) %>%
  mutate(label = factor(label, levels = taxa_order))

rec_linewidth <- case_when(
  n_taxa < 20  ~ 5,
  n_taxa < 50  ~ 2.5,
  n_taxa < 100 ~ 1.5,
  TRUE         ~ 0.8
)


###############################################################################
## 3. RECOMBINATION PANEL  (horizontal — original orientation)               ##
###############################################################################

p_rec <- ggplot(rec_long,
                aes(x = start, xend = end,
                    y = label, yend = label,
                    colour = rec_type)) +
  geom_segment(alpha = 0.6, linewidth = rec_linewidth) +
  scale_colour_manual(
    values = rec_colours,
    guide  = "none"          
  ) +
  scale_x_continuous(
    limits = c(region_start, region_end),
    expand = c(0, 0),
    labels = scales::comma
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Genomic position (bp)") +
  theme_bw() +
  theme(
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.y     = element_blank(),
    axis.line.y      = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(5, 10, 5, 0)
  )


###############################################################################
## 4. GENE ANNOTATION ARROWS (cmeABCR)                                       ##
###############################################################################

gene_plot_df <- gene_df %>%
  mutate(
    molecule = "cmeABCR_operon",
    midpoint = (start + end) / 2
  )

p_genes <- ggplot(gene_plot_df,
                  aes(xmin = start, xmax = end, y = molecule,
                      fill = gene, forward = direction == 1)) +
  geom_gene_arrow(
    arrowhead_height  = unit(4, "mm"),
    arrowhead_width   = unit(2, "mm"),
    arrow_body_height = unit(3, "mm")
  ) +
  geom_gene_label(aes(label = gene), align = "centre", grow = FALSE,
                  height = grid::unit(4, "mm"), fontface = "italic") +
  scale_fill_manual(values = gene_fills, guide = "none") +
  scale_x_continuous(limits = c(region_start, region_end), expand = c(0, 0)) +
  theme_genes() +
  theme(
    axis.text.y  = element_blank(),
    axis.title   = element_blank(),
    axis.ticks   = element_blank(),
    axis.line    = element_blank(),
    axis.text.x  = element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(2, 10, 0, 0)
  )


###############################################################################
## 5. RECOMBINATION DENSITY HEATMAP                                          ##
###############################################################################

density_df <- rec_region %>%
  mutate(bases = map2(start, end, seq)) %>%
  select(bases) %>%
  unnest(bases) %>%
  count(bases, name = "count") %>%
  mutate(
    count = pmin(count, quantile(count, 0.95)),
    y     = 1
  )

p_heatmap <- ggplot(density_df, aes(x = bases, y = y,
                                    fill = count, colour = count)) +
  geom_tile() +
  scale_x_continuous(limits = c(region_start, region_end), expand = c(0, 0)) +
  scale_fill_gradient2(
    low      = "navy", mid = "orange", high = "red",
    midpoint = max(density_df$count) / 2,
    name     = "Recombination\nevents",
    aesthetics = c("fill", "colour")
  ) +
  guides(
    fill   = guide_colorbar(title.position = "top", title.hjust = 0.5,
                            direction = "vertical", barheight = 4),
    colour = "none"
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "navy"),
    legend.position  = "right",
    plot.margin      = margin(0, 10, 0, 0)
  )


###############################################################################
## 6. METADATA PANEL                                                         ##
###############################################################################

meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Standardise id column name
if ("Id" %in% names(meta)) meta <- rename(meta, id = Id)
if ("ID" %in% names(meta)) meta <- rename(meta, id = ID)

meta_cols <- c("ST", "Source")  # <-- EDIT THESE
meta_cols <- meta_cols[meta_cols %in% names(meta)]

if (length(meta_cols) == 0) {
  warning("No matching metadata columns found. Skipping metadata panel.")
  p_meta <- NULL
  
} else {
  
  meta_long <- meta %>%
    select(id, all_of(meta_cols)) %>%
    filter(id %in% taxa_order) %>%
    mutate(across(all_of(meta_cols), as.character)) %>%
    pivot_longer(-id, names_to = "variable", values_to = "value") %>%
    mutate(
      id       = factor(id, levels = taxa_order),
      variable = factor(variable, levels = meta_cols),
      value    = as.character(value)
    )
  
  # ── Diagnostics ────────────────────────────────────────────────────────────
  missing_ST <- meta_long %>%
    filter(variable == "ST", is.na(value) | value == "" | value == "NA") %>%
    pull(id) %>%
    as.character()
  
  cat("Isolates with missing/blank ST:\n")
  if (length(missing_ST) == 0) cat("  (none)\n") else print(missing_ST)
  
  cat("\nAll ST values and counts:\n")
  meta_long %>%
    filter(variable == "ST") %>%
    count(value) %>%
    print(n = Inf)
  
  # ── Build colour palettes ───────────────────────────────────────────────────
  palette_names <- c("Set2", "Set3", "Paired", "Dark2", "Accent")
  meta_palettes <- list()
  
  for (i in seq_along(meta_cols)) {
    col    <- meta_cols[i]
    vals   <- unique(meta_long$value[meta_long$variable == col])
    n_vals <- length(vals)
    
    if (n_vals > 8) {
      colours <- many_colour_palette(n_vals)
    } else {
      pal_name <- palette_names[((i - 1) %% length(palette_names)) + 1]
      colours  <- brewer.pal(max(3, n_vals), pal_name)[seq_len(n_vals)]
    }
    
    names(colours)       <- vals
    meta_palettes[[col]] <- colours
  }
  
  # ── Build individual heatmap tile columns ───────────────────────────────────
  make_meta_tile <- function(col_name, pal) {
    df_sub <- meta_long %>% filter(variable == col_name)
    
    ggplot(df_sub, aes(x = variable, y = id, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_manual(values = pal, name = col_name, na.value = "grey80") +
      scale_y_discrete(drop = FALSE) +
      theme_void() +
      theme(
        axis.text.x     = element_text(angle = 45, hjust = 1, size = 8,
                                       face = "bold"),
        legend.position = "none",   # legends collected separately below
        plot.margin     = margin(5, 1, 5, 1)
      )
  }
  
  meta_plots <- map2(meta_cols, meta_palettes[meta_cols], make_meta_tile)
  p_meta     <- wrap_plots(meta_plots, nrow = 1)
}


###############################################################################
## 7. COMBINE ALL PANELS                                                     ##
###############################################################################

options("aplot_guides" = "keep")

# Main figure: tree + metadata + recombination, aligned on y (taxa)
if (!is.null(p_meta)) {
  combined <- p_rec %>%
    aplot::insert_left(p_meta, width = 0.06 * length(meta_cols)) %>%
    aplot::insert_left(p_tree, width = 0.4)
} else {
  combined <- p_rec %>%
    aplot::insert_left(p_tree, width = 0.4)
}

# Gene arrows and heatmap on top, aligned to recombination x-axis
combined <- combined %>%
  aplot::insert_top(p_genes,  height = 0.06) %>%
  aplot::insert_top(p_heatmap + theme(legend.position = "none"), height = 0.03)

# ── Collect legends ──────────────────────────────────────────────────────────

# Heatmap legend (vertical colour bar)
heatmap_legend <- cowplot::get_legend(
  p_heatmap + theme(legend.position = "right")
)

# Metadata legends (one per column, single-column layout)
if (!is.null(p_meta)) {
  meta_legends_explicit <- map2(meta_cols, meta_palettes[meta_cols], function(col_name, pal) {
    df_sub <- meta_long %>% filter(variable == col_name)
    
    p_leg <- ggplot(df_sub, aes(x = variable, y = id, fill = value)) +
      geom_tile() +
      scale_fill_manual(
        values   = pal,
        name     = col_name,
        na.value = "grey80",
        guide    = guide_legend(
          title.position = "top",
          title.hjust    = 0.5,
          ncol           = 1,
          override.aes   = list(size = 3)
        )
      ) +
      theme_void() +
      theme(
        legend.position = "right",
        legend.title    = element_text(face = "bold", size = 9),
        legend.text     = element_text(size = 8),
        legend.key.size = unit(0.4, "cm"),
        legend.key      = element_rect(colour = "white", linewidth = 0.3)
      )
    
    cowplot::get_legend(p_leg)
  })
  
  all_legends <- c(list(heatmap_legend), meta_legends_explicit)
} else {
  all_legends <- list(heatmap_legend)
}

all_legends <- Filter(Negate(is.null), all_legends)

# Stack all legends vertically on the right
legend_col <- cowplot::plot_grid(
  plotlist = all_legends,
  ncol     = 1,
  align    = "v",
  axis     = "lr"
)

# ── Final assembly ────────────────────────────────────────────────────────────
final_plot <- cowplot::plot_grid(
  aplot::as.patchwork(combined),
  legend_col,
  nrow       = 1,
  rel_widths = c(1, 0.2)   # increase 0.2 if legends are clipped
)

final_plot

###############################################################################
## 8. SAVE OUTPUT                                                            ##
###############################################################################

ggsave(
  filename = output_file,
  plot     = final_plot,
  width    = output_width,
  height   = output_height,
  dpi      = 600
)



###############################################################################
## NOTES                                                                     ##
## ─────────────────────────────────────────────────────────────────────────  ##
## 1. Gene coordinates:                                                      ##
##    - For C. coli: verify against your C. coli reference GFF               ##
##      grep -P "cme[ABCR]" your_reference.gff | awk '{print $4,$5,$7}'     ##
##                                                                           ##
## 2. Tip labels removed. Re-enable by uncommenting geom_tiplab() in        ##
##    Section 1 and increasing tree width back to 0.6 in Section 7.         ##
##                                                                           ##
## 3. If legends are clipped, increase rel_widths c(1, 0.2) → c(1, 0.3)    ##
##    in the final cowplot::plot_grid() call.                                ##
##                                                                           ##
## 4. Colour palettes auto-scale: ≤8 categories use RColorBrewer palettes;  ##
##    >8 categories use a 16-hue interpolated palette (many_colour_palette). ##
##                                                                           ##
## 5. Missing/blank ST values are printed to console and shown as grey80    ##
##    in the metadata tiles (na.value = "grey80").                           ##
##                                                                           ##
## 6. Metadata CSV format expected:                                          ##
##      id,ST,Source,Year                                                    ##
##      sample1,ST-48,Poultry,2018                                           ##
##      sample2,ST-45,Human,2019                                             ##
##                                                                           ##
## 7. The gene coordinates are for C. coli. For C. jejuni NCTC 11168:       ##
##      cmeR = Cj0368c, cmeA = Cj0367c, cmeB = Cj0366c, cmeC = Cj0365c    ##
###############################################################################
