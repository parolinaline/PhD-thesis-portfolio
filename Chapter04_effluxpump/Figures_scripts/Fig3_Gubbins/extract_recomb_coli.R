

getwd()
setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 04 Campy cmeABC efflux pump/WRITING/Figures/FIG_gubbins/")

###############################################################################
## Extract recombination events overlapping cmeABCR genes                    ##
## Output: CSV table with isolate, recombination coords, and affected genes  ##
###############################################################################

library(tidyverse)

# ── INPUT FILES ──────────────────────────────────────────────────────────────
rec_gff_file <- "coli_gubbins_recombination.gff"
output_csv   <- "cmeABCRcoli_recombination_table.csv"

# ── GENE COORDINATES (C. coli) ───────────────────────────────────────────────
gene_df <- tibble::tribble(
  ~gene,   ~gene_start, ~gene_end,
  "cmeC",  373168,      374646,
  "cmeB",  374639,      377761,
  "cmeA",  377761,      378864,
  "cmeR",  378972,      379607
)

# ── PARSE GUBBINS GFF ───────────────────────────────────────────────────────
raw <- readLines(rec_gff_file)
raw <- raw[!grepl("^##", raw)]

rec_df <- read.table(text = raw, sep = "\t", header = FALSE,
                     stringsAsFactors = FALSE, fill = TRUE, quote = "") %>%
  setNames(c("seq", "prog", "type", "start", "end",
             "score", "strand", "frame", "attributes")) %>%
  rowwise() %>%
  mutate(
    taxa_str = str_extract(attributes, "taxa=\\s*[^;]+"),
    taxa_str = str_remove(taxa_str, "taxa=\\s*"),
    taxa_str = str_trim(taxa_str),
    # FIX: strip all quote characters, split, and remove empty strings
    taxa     = list(
      str_split(taxa_str, "\\s+")[[1]] %>%
        str_remove_all('"') %>%
        .[nchar(.) > 0]
    ),
    n_taxa   = length(taxa),
    rec_type = if_else(n_taxa == 1, "unique", "shared")
  ) %>%
  ungroup() %>%
  select(start, end, taxa, n_taxa, rec_type)

# ── FIND OVERLAPS WITH cmeABCR GENES ────────────────────────────────────────
region_start <- min(gene_df$gene_start)
region_end   <- max(gene_df$gene_end)

rec_in_region <- rec_df %>%
  filter(end >= region_start & start <= region_end)

rec_with_genes <- rec_in_region %>%
  rowwise() %>%
  mutate(
    affected_genes = list(
      gene_df %>%
        filter(gene_end >= start & gene_start <= end) %>%
        pull(gene)
    ),
    genes_hit = paste(affected_genes, collapse = ", "),
    n_genes   = length(affected_genes)
  ) %>%
  ungroup()

# Unnest to one row per taxon
result_table <- rec_with_genes %>%
  unnest(taxa) %>%
  select(
    isolate       = taxa,
    rec_start     = start,
    rec_end       = end,
    rec_type,
    n_taxa_in_rec = n_taxa,
    genes_hit,
    n_genes
  ) %>%
  mutate(rec_length = rec_end - rec_start + 1) %>%
  arrange(isolate, rec_start)

# ── SANITY CHECK: confirm no stray quotes remain ─────────────────────────────
stray <- result_table %>% filter(str_detect(isolate, '"'))
if (nrow(stray) > 0) {
  warning("Stray quote characters still detected in isolate names!")
  print(stray)
} else {
  cat("✓ No stray quotes detected in isolate names.\n")
}

# ── SUMMARY TABLE: one row per isolate ──────────────────────────────────────
summary_table <- result_table %>%
  group_by(isolate) %>%
  summarise(
    n_recombinations = n(),
    total_rec_bp     = sum(rec_length),
    has_cmeC         = any(str_detect(genes_hit, "cmeC")),
    has_cmeB         = any(str_detect(genes_hit, "cmeB")),
    has_cmeA         = any(str_detect(genes_hit, "cmeA")),
    has_cmeR         = any(str_detect(genes_hit, "cmeR")),
    all_genes_hit    = paste(unique(unlist(str_split(genes_hit, ", "))),
                             collapse = ", "),
    rec_types        = paste(unique(rec_type), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_recombinations))

# ── SAVE OUTPUTS ────────────────────────────────────────────────────────────
write.csv(result_table,
          file = output_csv,
          row.names = FALSE)

write.csv(summary_table,
          file = gsub(".csv", "_summary.csv", output_csv),
          row.names = FALSE)

# ── PRINT SUMMARY ───────────────────────────────────────────────────────────
cat("\n=== cmeABCR Recombination Summary ===\n\n")
cat("Total recombination events in region:", nrow(rec_in_region), "\n")
cat("Isolates affected:", n_distinct(result_table$isolate), "\n\n")

cat("Per gene:\n")
cat("  cmeC:", sum(summary_table$has_cmeC), "isolates\n")
cat("  cmeB:", sum(summary_table$has_cmeB), "isolates\n")
cat("  cmeA:", sum(summary_table$has_cmeA), "isolates\n")
cat("  cmeR:", sum(summary_table$has_cmeR), "isolates\n\n")

