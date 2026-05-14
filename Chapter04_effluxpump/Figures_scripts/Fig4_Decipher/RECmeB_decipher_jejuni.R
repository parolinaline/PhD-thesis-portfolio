#!/usr/bin/env Rscript
# RE_cmeB_alignment_visualization.R
# Visualize RE-CmeB mutations using DECIPHER alignment
# For thesis Chapter 4 figures

library(DECIPHER)
library(Biostrings)
library(ggplot2)
library(dplyr)

setwd("C:/Users/21004751/OneDrive - Massey University/Student Projects - Aline Parolin Calarga PhD Project - Aline Parolin Calarga PhD Project/Thesis/Chapter 04 Campy cmeABC efflux pump/WRITING/Figures/RE_CmeB/")

input_file <- "RE_cmeB_jej_7aa.faa"

cat("RE-CmeB Alignment Visualization\n")
cat("===============================\n\n")

cat("Loading RE-CmeB sequences from:", input_file, "\n")
aa_sequences <- readAAStringSet(input_file)

cat("Loaded", length(aa_sequences), "amino acid sequences\n")
cat("Sequence names:\n")
print(names(aa_sequences))

cat("\nAligning sequences using DECIPHER...\n")
aligned_sequences <- AlignSeqs(aa_sequences, verbose = TRUE)

cat("Alignment complete!\n")
cat("Alignment length:", width(aligned_sequences)[1], "amino acids\n\n")

# ============================================================
# All positions to track (expanded for new combinations)
# ============================================================
TARGET_POSITIONS <- c(136, 140, 257, 267, 292, 294, 303)

# Expected mutant amino acids at each position
MUTANT_AA <- c(
  "136" = "A",   # T136A
  "140" = "G",   # A140G
  "257" = "D",   # N257D
  "267" = "I",   # V267I
  "292" = "I",   # M292I
  "294" = "D",   # N294D
  "303" = "N"    # H303N
)

# ============================================================
# Annotate mutation status
# ============================================================
annotate_mutation_status <- function(aligned_seqs) {
  
  # Initialise annotation data frame
  pos_cols <- paste0("pos_", TARGET_POSITIONS)
  annotations <- data.frame(
    sequence_id    = names(aligned_seqs),
    matrix(NA_character_, nrow = length(aligned_seqs), ncol = length(TARGET_POSITIONS),
           dimnames = list(NULL, pos_cols)),
    mutation_status = NA_character_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(aligned_seqs)) {
    ungapped_seq <- gsub("-", "", as.character(aligned_seqs[[i]]))
    seq_len      <- nchar(ungapped_seq)
    
    # Extract amino acid at each target position
    for (pos in TARGET_POSITIONS) {
      col <- paste0("pos_", pos)
      if (seq_len >= pos) {
        annotations[[col]][i] <- substr(ungapped_seq, pos, pos)
      }
    }
    
    # Helper: is a given mutation present?
    has_mut <- function(pos) {
      col <- paste0("pos_", pos)
      aa  <- annotations[[col]][i]
      !is.na(aa) && aa == MUTANT_AA[as.character(pos)]
    }
    
    # Classify — most-complete combination first
    annotations$mutation_status[i] <- dplyr::case_when(
      # Full 7-mutation RE-CmeB
      has_mut(136) & has_mut(140) & has_mut(257) & has_mut(267) &
        has_mut(292) & has_mut(294) & has_mut(303) ~
        "T136A+A140G+N257D+V267I+M292I+N294D+H303N",
      
      # 4-mutation combination
      has_mut(136) & has_mut(257) & has_mut(267) & has_mut(292) ~
        "T136A+N257D+V267I+M292I",
      
      # Original complete RE-CmeB (2-mutation)
      has_mut(136) & has_mut(292) ~
        "T136A+M292I",
      
      # Legacy partial combinations
      has_mut(136) ~ "T136A only",
      has_mut(292) ~ "M292I only",
      
      TRUE ~ "Wild-type or other variant"
    )
  }
  
  return(annotations)
}

annotations <- annotate_mutation_status(aligned_sequences)

# Print summary
cat("=== MUTATION SUMMARY ===\n")
print(table(annotations$mutation_status))

cat("\nPer-position variant counts:\n")
for (pos in TARGET_POSITIONS) {
  col <- paste0("pos_", pos)
  cat(sprintf("\nPosition %d:\n", pos))
  print(table(annotations[[col]], useNA = "ifany"))
}

# ============================================================
# Find alignment positions for all target sites
# ============================================================
find_alignment_positions <- function(aligned_seqs, target_positions = TARGET_POSITIONS) {
  ref_seq <- as.character(aligned_seqs[[1]])
  position_map  <- list()
  ungapped_pos  <- 0L
  
  for (align_pos in seq_len(nchar(ref_seq))) {
    if (substr(ref_seq, align_pos, align_pos) != "-") {
      ungapped_pos <- ungapped_pos + 1L
      if (ungapped_pos %in% target_positions) {
        position_map[[as.character(ungapped_pos)]] <- align_pos
      }
    }
  }
  return(position_map)
}

alignment_positions <- find_alignment_positions(aligned_sequences)

cat("\n=== ALIGNMENT POSITION MAPPING ===\n")
for (pos in names(alignment_positions)) {
  cat(sprintf("Original position %s → alignment position %d\n",
              pos, alignment_positions[[pos]]))
}

# ============================================================
# DECIPHER HTML outputs
# ============================================================
cat("\n=== CREATING ALIGNMENT VISUALIZATIONS (HTML) ===\n")

full_html <- "RE_cmeB_full_alignment.html"
BrowseSeqs(aligned_sequences, htmlFile = full_html, openURL = FALSE)
cat("Full alignment saved:", full_html, "\n")

if (!is.null(alignment_positions[["136"]]) && !is.null(alignment_positions[["303"]])) {
  
  start_pos  <- alignment_positions[["136"]]
  end_pos    <- alignment_positions[["303"]]
  
  region_seqs <- subseq(aligned_sequences, start = start_pos, end = end_pos)
  
  # Rename sequences to make the position range explicit in the viewer
  names(region_seqs) <- paste0(names(region_seqs), "  [CmeB positions 136-303]")
  
  region_html <- "RE_cmeB_region_136_303.html"
  
  BrowseSeqs(region_seqs,
             htmlFile = region_html,
             openURL  = FALSE,
             title    = "CmeB alignment — positions 136 to 303 (RE-CmeB mutation sites)")
  
  cat("Focused alignment (positions 136–303) saved:", region_html, "\n")
  
} else {
  cat("WARNING: Could not map positions 136 and/or 303 from the reference sequence.",
      "Check that your sequences are long enough.\n")
}

# ============================================================
# Summary bar chart (ggplot2)
# ============================================================
cat("\n=== CREATING MUTATION SUMMARY FIGURE (PDF) ===\n")

mutation_counts <- as.data.frame(table(annotations$mutation_status)) %>%
  rename(mutation_status = Var1, count = Freq)

status_colours <- c(
  "T136A+A140G+N257D+V267I+M292I+N294D+H303N" = "#7b2d8b",   # purple  — 7-mut
  "T136A+N257D+V267I+M292I"                    = "#d62728",   # red     — 4-mut
  "T136A+M292I"                                = "#ff7f0e",   # orange  — 2-mut
  "T136A only"                                 = "#ffbb78",   # light orange
  "M292I only"                                 = "#aec7e8",   # light blue
  "Wild-type or other variant"                 = "#1f77b4"    # blue
)

p <- ggplot(mutation_counts,
            aes(x = reorder(mutation_status, count), y = count, fill = mutation_status)) +
  geom_col(width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = status_colours) +
  labs(title = "RE-CmeB Mutation Status", x = NULL, y = "Number of sequences") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("RE_cmeB_mutation_TEST.pdf", p, width = 10, height = 5)
cat("Mutation summary figure saved: RE_cmeB_mutation_summary.pdf\n")

# ============================================================
# Mutation table
# ============================================================
create_mutation_table <- function(annotations) {
  annotations %>%
    mutate(
      Resistance_Enhancement = dplyr::case_when(
        mutation_status == "T136A+A140G+N257D+V267I+M292I+N294D+H303N" ~ "Very High",
        mutation_status == "T136A+N257D+V267I+M292I"                   ~ "High",
        mutation_status == "T136A+M292I"                               ~ "Moderate",
        grepl("only", mutation_status)                                 ~ "Low",
        TRUE                                                           ~ "None"
      )
    ) %>%
    select(sequence_id, starts_with("pos_"), mutation_status, Resistance_Enhancement)
}

mutation_table <- create_mutation_table(annotations)

# ============================================================
# Save all outputs
# ============================================================
cat("\n=== SAVING TABULAR RESULTS ===\n")

writeXStringSet(aligned_sequences, "RE_cmeB_aligned_sequencesNEW.fasta")
cat("Aligned FASTA saved: RE_cmeB_aligned_sequences.fasta\n")

write.csv(mutation_table,  "RE_cmeB_mutation_tableNEW.csv",   row.names = FALSE)
cat("Mutation table saved: RE_cmeB_mutation_table.csv\n")

write.csv(annotations, "RE_cmeB_annotationsNEW.csv", row.names = FALSE)
cat("Full annotations saved: RE_cmeB_annotations.csv\n")

cat("\n=== ALL DONE ===\n")
cat("HTML alignment files can be opened in any web browser.\n")
