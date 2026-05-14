
getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 04 Campy cmeABC efflux pump/WRITING/Figures/FIG_gubbins/")

###############################################################################
## Visualise Gubbins recombination in the cmeABCR region                     ##
## Aligned to phylogenetic tree with gene annotation, metadata & heatmap     ##
##                                                                           ##
## Usage:                                                                    ##
##   Rscript plot_gubbins_cmeABCR.R                                          ##
##                                                                           ##
## IMPORTANT: Before running, update the following:                          ##
##   1. File paths in the "INPUT FILES" section                              ##
##   2. Gene coordinates in the "GENE COORDINATES" section                   ##
##      (extract from your reference GFF with grep "cme" ref.gff)            ##
##   3. Metadata column names in the "METADATA" section                      ##
##   4. Output file names                                                    ##
###############################################################################

# ── Libraries ────────────────────────────────────────────────────────────────
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(gggenes)      # for gene arrow diagrams
library(aplot)
library(patchwork)
library(cowplot)
library(RColorBrewer)


# ── INPUT FILES ──────────────────────────────────────────────────────────────
tree_file     <- "./gubbins_jejuni_tree.tre"
rec_gff_file  <- "gubbins_jejuni_predictions.gff"
meta_file     <- "Metadata_jejuni.csv"   # CSV with columns: id, ST, Source, Year (etc.)
output_file   <- "jejuni_cmeABCR_recombination.pdf"


# ── GENE COORDINATES ────────────────────────────────────────────────────────
# *** VERIFY THESE FROM YOUR REFERENCE GFF ***
#   grep -i "cme" your_reference.gff
#
# Using reference genome SC0017 hybrid assembly
# C. jejuni
gene_df <- tibble::tribble(
  ~gene,   ~start,  ~end,    ~strand, ~direction,
  "cmeC",  342811,  344289,  -1,      -1,
  "cmeB",  344282,  347404,  -1,      -1,
  "cmeA",  347404,  348507,  -1,      -1,
  "cmeR",  348601,  349233,  -1,      -1
)

# Define the genomic window (add some padding around the operon)
padding       <- 500
region_start  <- min(gene_df$start) - padding
region_end    <- max(gene_df$end)   + padding

# ── OUTPUT SETTINGS ──────────────────────────────────────────────────────────
output_width  <- 14    # inches
output_height <- 10    # inches

# ── COLOUR PALETTES ─────────────────────────────────────────────────────────
rec_colours <- c("shared" = "#D73027", "unique" = "#4575B4")
gene_fills  <- c("cmeR" = "#FDB462", "cmeA" = "#80B1D3",
                 "cmeB" = "#FB8072", "cmeC" = "#B3DE69")

# ── HELPER: scalable palette generator ───────────────────────────────────────
# Returns n colours by interpolating a named RColorBrewer palette.
# Works for any number of categories (including >12).
make_palette <- function(n, brewer_pal = "Spectral") {
  base_cols <- brewer.pal(brewer.pal.info[brewer_pal, "maxcolors"], brewer_pal)
  colorRampPalette(base_cols)(n)
}

###############################################################################
## 1. READ AND PROCESS TREE                                                  ##
###############################################################################

tree   <- ape::read.tree(tree_file)
n_taxa <- length(tree$tip.label)

# Adjust linewidth based on number of taxa
tree_linewidth <- case_when(
  n_taxa < 20  ~ 0.8,
  n_taxa < 50  ~ 0.5,
  n_taxa < 100 ~ 0.35,
  TRUE          ~ 0.25
)

p_tree <- ggtree(tree, linewidth = tree_linewidth) +   # linewidth replaces deprecated size
  theme_tree2() +
  theme(
    legend.position = "none",
    axis.line.x     = element_line(linewidth = 0.25),
    plot.margin     = margin(5, 0, 5, 5)
  )

# Get the taxa order from the tree (important for aligning all panels)
taxa_order <- rev(get_taxa_name(p_tree))

p_tree

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

rec_df <- parse_gubbins_gff(rec_gff_file)

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
  mutate(label = factor(label, levels = taxa_order))

# Adjust segment linewidth based on number of taxa
rec_linewidth <- case_when(
  n_taxa < 20  ~ 5,
  n_taxa < 50  ~ 2.5,
  n_taxa < 100 ~ 1.5,
  TRUE          ~ 0.8
)

###############################################################################
## 3. RECOMBINATION PANEL                                                    ##
###############################################################################

p_rec <- ggplot(rec_long,
                aes(x = start, xend = end,
                    y = label, yend = label,
                    colour = rec_type)) +
  geom_segment(alpha = 0.6, linewidth = rec_linewidth) +
  scale_colour_manual(
    values = rec_colours,
    name   = "Recombination\ntype"
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
    legend.position  = "right",
    plot.margin      = margin(5, 10, 5, 0)
  )

###############################################################################
## 4. GENE ANNOTATION ARROWS (cmeABCR)                                       ##
###############################################################################

gene_plot_df <- gene_df %>%
  mutate(
    molecule = "cmeABCR_operon",
    xmin     = if_else(direction == -1, end, start),
    xmax     = if_else(direction == -1, start, end),
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
  scale_x_continuous(
    limits = c(region_start, region_end),
    expand = c(0, 0)
  ) +
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
    count = pmin(count, quantile(count, 0.95)),  # cap at 95th percentile
    y = 1
  )

p_heatmap <- ggplot(density_df, aes(x = bases, y = y,
                                    fill = count, colour = count)) +
  geom_tile() +
  scale_x_continuous(
    limits = c(region_start, region_end),
    expand = c(0, 0)
  ) +
  scale_fill_gradient2(
    low      = "navy",
    mid      = "orange",
    high     = "red",
    midpoint = max(density_df$count) / 2,
    name     = "Recombination\nevents",
    aesthetics = c("fill", "colour")
  ) +
  guides(
    fill   = guide_colorbar(title.position = "top",
                            title.hjust    = 0.5,
                            direction      = "vertical",
                            barheight      = 6),
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

# *** UPDATE: select the metadata columns you want to display ***
meta_cols <- c("ST", "Source")  # <-- EDIT THESE

# Keep only columns that exist in the file
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
  
  # ── Colour palettes: one per metadata column ────────────────────────────
  # Each column gets a different base palette so ST and Source are visually
  # distinct. colorRampPalette interpolates to however many categories exist,
  # so it works regardless of how many STs you have.
  base_palettes <- c(
    "Spectral",   # ST  – wide spectral range, great for many categories
    "Set2",       # Source – qualitative, good for a handful of sources
    "Paired",     # extra columns if needed
    "Dark2",
    "Set3"
  )
  
  meta_palettes <- list()
  
  for (i in seq_along(meta_cols)) {
    col      <- meta_cols[i]
    vals     <- sort(unique(meta_long$value[meta_long$variable == col]))
    n_vals   <- length(vals)
    pal_name <- base_palettes[((i - 1) %% length(base_palettes)) + 1]
    colours  <- make_palette(n_vals, pal_name)
    names(colours) <- vals
    meta_palettes[[col]] <- colours
  }
  
  # ── Build one tile panel per metadata column ────────────────────────────
  make_meta_tile <- function(col_name, pal) {
    df_sub <- meta_long %>% filter(variable == col_name)
    
    ggplot(df_sub, aes(x = col_name, y = id, fill = value)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      scale_fill_manual(values = pal, name = col_name, na.value = "grey80") +
      scale_y_discrete(drop = FALSE) +
      theme_void() +
      theme(
        axis.text.x    = element_text(angle = 45, hjust = 1, size = 8,
                                      face = "bold"),
        legend.position = "right",
        legend.title    = element_text(face = "bold", size = 8),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.3, "cm"),
        plot.margin     = margin(5, 1, 5, 1)
      ) +
      guides(fill = guide_legend(
        title.position = "top",
        ncol           = min(6, length(pal)),   # up to 6 columns in legend
        title.hjust    = 0.5
      ))
  }
  
  meta_plots <- map2(meta_cols, meta_palettes[meta_cols], make_meta_tile)
  
  # Strip legends from panels (collected separately below)
  meta_plots_no_legend <- map(meta_plots, ~ .x + theme(legend.position = "none"))
  meta_legends         <- map(meta_plots, ~ cowplot::get_legend(.x))
  
  p_meta <- wrap_plots(meta_plots_no_legend, nrow = 1)
}

###############################################################################
## 7. COMBINE ALL PANELS                                                     ##
###############################################################################

options("aplot_guides" = "keep")

if (!is.null(p_meta)) {
  combined <- p_rec %>%
    aplot::insert_left(p_meta, width = 0.15 * length(meta_cols)) %>%
    aplot::insert_left(p_tree, width = 0.4)
} else {
  combined <- p_rec %>%
    aplot::insert_left(p_tree, width = 0.4)
}

combined <- combined %>%
  aplot::insert_top(p_genes,  height = 0.06) %>%
  aplot::insert_top(p_heatmap + theme(legend.position = "none"), height = 0.03)

# ── Collect legends ──────────────────────────────────────────────────────────
rec_legend     <- cowplot::get_legend(p_rec)
heatmap_legend <- cowplot::get_legend(p_heatmap)

all_legends <- if (!is.null(p_meta)) {
  c(list(rec_legend, heatmap_legend), meta_legends)
} else {
  list(rec_legend, heatmap_legend)
}

legend_col <- cowplot::plot_grid(plotlist = all_legends, ncol = 1,
                                 align = "v", axis = "l")

# ── Final assembly ────────────────────────────────────────────────────────────
final_plot <- cowplot::plot_grid(
  aplot::as.patchwork(combined),
  legend_col,
  nrow        = 1,
  rel_widths  = c(1, 0.18)
)

final_plot

head(tree$tip.label, 10)

# Full diagnostic
cat("Tree tips not in metadata:\n")
print(setdiff(taxa_order, meta$id))

cat("\nMetadata IDs not in tree:\n")
print(setdiff(meta$id, taxa_order))

cat("\nRows with NA in ST or Source after join:\n")
meta %>% filter(id %in% taxa_order) %>% filter(is.na(ST) | is.na(Source)) %>% print()

###############################################################################
## 8. SAVE OUTPUT                                                            ##
###############################################################################

ggsave(
  filename = output_file,
  plot     = final_plot,
  width    = output_width,
  height   = output_height,
  dpi      = 300
)

cat("✓ Plot saved to:", output_file, "\n")
cat("  Region:", scales::comma(region_start), "-",
    scales::comma(region_end), "bp\n")
cat("  Taxa:", n_taxa, "\n")
cat("  Recombination events in region:", nrow(rec_region), "\n")

###############################################################################
## NOTES                                                                     ##
## ─────────────────────────────────────────────────────────────────────────  ##
## 1. Gene coordinates:                                                      ##
##    - For C. jejuni: verify against your reference GFF                     ##
##    - For C. coli: extract from your C. coli reference                     ##
##      grep -P "cme[ABCR]" your_reference.gff | awk '{print $4,$5,$7}'     ##
##                                                                           ##
## 2. If gggenes is not installed:                                           ##
##      install.packages("gggenes")                                          ##
##                                                                           ##
## 3. Metadata CSV format expected:                                          ##
##      id,ST,Source,Year                                                    ##
##      sample1,ST-48,Poultry,2018                                           ##
##      sample2,ST-45,Human,2019                                             ##
##                                                                           ##
## 4. To compare C. jejuni vs C. coli side by side, run this script twice    ##
##    with the respective inputs, then combine PDFs or use cowplot::plot_grid##
##    on the two final plots.                                                ##
##                                                                           ##
## 5. The gene coordinates below are for NCTC 11168 (AL111168).              ##
##    In the cmeABC operon, genes are on the complement strand:              ##
##      cmeR = Cj0368c                                                       ##
##      cmeA = Cj0367c (periplasmic fusion protein)                          ##
##      cmeB = Cj0366c (inner membrane RND transporter)                      ##
##      cmeC = Cj0365c (outer membrane channel)                              ##
##    Note: cmeR is upstream and transcribed divergently from cmeABC.        ##
###############################################################################